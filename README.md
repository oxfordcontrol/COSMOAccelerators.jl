# COSMOAccelerators.jl

## Description
`COSMOAccelerators` defines an abstract type `AbstractAccelerator` that can be used as an interface to write accelerator methods for algorithms based on non-expansive operators, e.g. fixed-point iterations. 
The acceleration methods in this package were originally developed for the convex conic solver [COSMO.jl](https://github.com/oxfordcontrol/COSMO.jl), but can also be used with other fixed-point methods, e.g. algorithms based on the [ProximalOperators.jl](https://github.com/kul-forbes/ProximalOperators.jl) package.

## Installation

`COSMOAccelerators` can be added via the Julia package manager (type `]`): `pkg> add COSMOAccelerators`

## Usage
The `AbstractAccelerator` interface assumes that the solver method `T` iteratively updates a vector `g = T(x)`. The interface for an abstract accelerator `aa` requires implementation of the following methods:

| Methods to implement        | Brief description           |  
| ------------- |:-------------| 
| `update!(aa, g, x, iter)`      | Use iteration pairs to update the accelerator history | 
| `accelerate!(g, x, aa, iter)`     | Recombine past iterates and overwrite `g` with the accelerated point      |  
| `was_successful!(aa)` |  Tell the solver algorithms whether the last acceleration attempt was successful, e.g. no numerical issues   |
| `restart!(aa)` |        Restart the accelerator, i.e. forget any state or history |


| Optional methods        | Default definition | Brief description           |  
| ------------- |:-------------| :-------------|
| `log!(aa, args...; kwargs...)` |  do nothing     |  Algorithm tells accelerator to log something |

## Accelerators
The following accelerator types are currently exported:

- `EmptyAccelerator`: An accelerator that does nothing
- `AndersonAccelerator{T, BT, MT, RE}(dim; mem)`: implements a variant of a limited-memory Anderson Acceleration method
  - `T <: AbstractFloat`
  - `BT <: AbstractBroydenType`:
    - `Type1` for Anderson-Acceleration with a Broyden type1-update
    - `Type2{NormalEquations}` for Anderson-Acceleration with a Broyden type2-update and where the least squares subproblem is solved with normal equations
    - `Type2{QRDecomp}` for Anderson-Acceleration with a Broyden type2-update and where the least squares subproblem is solved with a updated QR decomposition
  - `MT <: AbstractMemory`:
    - `RollingMemory`: Once the memory limit of the history is reached, new iterates overwrite the oldest iterates
    - `RestartedMemory`: Once the memory limit of the history is reached, the history is cleared and the method starts from scratch
  - `RE <: AbstractRegularizer`:
    - `TikonovRegularizer`: Tikonov regularisation of the least-squares subproblem to improve conditioning 
    - `NoRegularizer`: No regularisation of the least-squares subproblem 

## Example
This is a modifed LASSO example from [ProximalAlgorithms.jl](https://github.com/kul-forbes/ProximalAlgorithms.jl/blob/master/test/problems/test_lasso_small.jl).
The first function defines the Douglas-Rachford algorithm based on proximal operators.
The second function adds an Anderson Acceleration step into the algorithm loop by passing the iterates to `update!` and `accelerate!`. The number of iterations for convergence decreases from 36 to 12.
```julia
using COSMOAccelerators 
using ProximalOperators, ProximalAlgorithms
using LinearAlgebra, Test

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

# define functions for proximal operators
f = LeastSquares(A, b)
g = NormL1(lam)

# optimal solution 
x_star = T[-3.877278911564627e-01, 0, 0, 2.174149659863943e-02, 6.168435374149660e-01]

# Define Douglas-Rachford splitting
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

# Define Douglas-Rachford splitting and Anderson Acceleration
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



```

## Related packages
- [COSMO.jl](https://github.com/oxfordcontrol/COSMO.jl) - an ADMM based solver for convex conic optimisation problems.
- [ProximalOperators.jl](https://github.com/kul-forbes/ProximalOperators.jl) - a set of operator building blocks for proximal algorithms

## Credits
`COSMOAccelerators.jl` is developed by [Michael Garstka](www.migarstka.com) at the University of Oxford.

