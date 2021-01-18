using COSMOAccelerators 
using ProximalOperators
using ProximalAlgorithms
using LinearAlgebra
using Test

# problem taken from 
# https://github.com/kul-forbes/ProximalAlgorithms.jl/blob/master/test/problems/test_lasso_small.jl
T = Float64
A = T[
    1.0 -2.0 3.0 -4.0 5.0
    2.0 -1.0 0.0 -1.0 3.0
    -1.0 0.0 4.0 -3.0 2.0
    -1.0 -1.0 -1.0 1.0 3.0
]

b = T[1.0, 2.0, 3.0, 4.0]
m, n = size(A)
R = real(T)
lam = R(0.1) * norm(A' * b, Inf)

# functions for proximal operator
f = LeastSquares(A, b)
g = NormL1(lam)

# correct solution 
x_star = T[-3.877278911564627e-01, 0, 0, 2.174149659863943e-02, 6.168435374149660e-01]

# Douglas-Rachford splitting
function solve(iter::ProximalAlgorithms.DRS_iterable, maxiter, gamma, tol)
    tol_stop(state::ProximalAlgorithms.DRS_state) = norm(state.res, Inf) / gamma <= tol

    iter = ProximalAlgorithms.halt(iter, tol_stop) # tolerance check
    iter = Base.Iterators.take(iter, maxiter) # max iter
    iter = enumerate(iter) # state = (k, state)
    
    state_final = nothing
    for state in iter
        state_final = state
    end
    return state_final[2].y, state_final[2].x, state_final[1]
end

# Douglas-Rachford splitting and Anderson Acceleration
function solve_and_accelerate(iter::ProximalAlgorithms.DRS_iterable, aa::AndersonAccelerator, maxiter, gamma, tol)
    tol_stop(state::ProximalAlgorithms.DRS_state) = norm(state.res, Inf) / gamma <= tol
    x_prev = copy(iter.x)

    iter = ProximalAlgorithms.halt(iter, tol_stop)
    iter = Base.Iterators.take(iter, maxiter) 
    
    state_final = nothing
    k_final = nothing
    for (k, state) in enumerate(iter)
        state_final = state
        k_final = k
    
        # accelerate
        update!(aa, state.x, x_prev, k)
        accelerate!(state.x, x_prev, aa, k)
        @. x_prev = state.x

    end
    return state_final.y, state_final.z, k_final
end



# test
@testset "Accelerator and Douglas-Rachford (ProximalOperators.jl)" begin
    tol = 1e-6
    maxiter = 200
    x0 = zeros(5)
    gamma = R(10.0) / opnorm(A)^2

    # vanilla Douglas-Rachford
    iter = ProximalAlgorithms.DRS_iterable(f, g, x0, gamma)
    y, z, it = solve(iter, maxiter, gamma, tol)
    @test norm(y - x_star, Inf) <= tol


    # Douglas-Rachford + Anderson Acceleration
    iter = ProximalAlgorithms.DRS_iterable(f, g, x0, gamma)
    aa = AndersonAccelerator(length(x0))
    y_aa, z_aa, it_aa = solve_and_accelerate(iter, aa, maxiter, gamma, tol)
    @test norm(y_aa - x_star, Inf) <= tol
    @test it_aa < it 

end


