export NoRegularizer, TikonovRegularizer, Type1, Type2, RollingMemory, RestartedMemory, AndersonAccelerator
export QRDecomp, NormalEquations

# ---------------------------
# AndersonAccelerator
# ---------------------------

"""
    AbstractRegularizer

Abstract supertype for Anderson Acceleration Type-II regularisation schemes.
"""
abstract type AbstractRegularizer end
struct TikonovRegularizer <: AbstractRegularizer end
# struct FrobeniusNormRegularizer <: AbstractRegularizer end
struct NoRegularizer <: AbstractRegularizer end

"""
    AbstractLeastSquaresMethod

Abstract supertype for least squares solution method when Type-II acceleration is used
"""
abstract type AbstractLeastSquaresMethod end
struct NormalEquations <: AbstractLeastSquaresMethod end
struct QRDecomp <: AbstractLeastSquaresMethod end

"""
    AbstractBroydenType

Abstract supertype for the Broyden type of the accelerator.
"""
abstract type AbstractBroydenType end
struct Type1 <: AbstractBroydenType end
struct Type2{LSM} <: AbstractBroydenType 
  function Type2{LSM}() where {LSM <: AbstractLeastSquaresMethod}
    return new{LSM}()
  end
end
Type2() = Type2{QRDecomp}()

"""
    AbstractMemory

Abstract supertype for the memory management of the accelerator.
"""
abstract type AbstractMemory end

"""
    RollingMemory

The accelerator will append new iterates to the history and discard the oldest iterate.
"""
struct RollingMemory <: AbstractMemory end

"""
    RestartedMemory

The accelerator will delete the history once the memory buffers are full.
"""
struct RestartedMemory <: AbstractMemory end

"""
    AndersonAccelerator{T, BT, MT, RE} <: AbstractAccelerator

Accelerator object implementing Anderson Acceleration. Parameterized by:

 - T: AbstractFloat, floating-point type
 - BT: Broyden-type, i.e. Type1 or Type2
 - MT: AbstractMemory, how full memory buffers are handled
 - RE: AbstractRegularizer

"""
mutable struct AndersonAccelerator{T, BT, MT, RE}  <: AbstractAccelerator
  init_phase::Bool
  mem::Int
  min_mem::Int
  dim::Int
  iter::Int
  num_accelerated_steps::Int
  fail_counter::Array{Int}
  fail_eta::Array{Int}
  fail_singular::Array{Int}
  x_last::Vector{T}
  g_last::Vector{T}
  f::Vector{T}
  f_last::Vector{T}
  eta::Vector{T}
  F::Matrix{T}
  X::Matrix{T}
  G::Matrix{T}
  M::Matrix{T}
  Q::Matrix{T}
  R::Matrix{T}
  λ::T # regularisation parameter
  activated::Bool # a flag that shows that the accelerator has been started
  success::Bool # a flag to indicate whether the last attempted acceleration was successful
  update_time::Float64
  accelerate_time::Float64
  acceleration_status::Vector{Tuple{Int, Symbol}}
  safeguarding_status::Vector{Tuple{Int, T, T, T}} 
  num_safe_accepted::Int
  num_safe_declined::Int
  activate_logging::Bool
  
  function AndersonAccelerator{T, BT, MT, RE}() where {T <: AbstractFloat, RE <: AbstractRegularizer, BT <: AbstractBroydenType, MT <: AbstractMemory}
    new(true, 0, 3, 0, 0, 0, zeros(Int, 0), zeros(Int, 0), zeros(Int, 0), zeros(T, 1), zeros(T, 1), zeros(T, 1),  zeros(T, 1), zeros(T, 1), zeros(T, 1), zeros(T, 1, 1), zeros(T, 1, 1), zeros(T, 1, 1), zeros(T, 1, 1), zeros(T, 1, 1), zeros(T, 1, 1), zero(T), false, false, 0., 0., Vector{Tuple{Int, Symbol}}(undef, 0),  Vector{Tuple{Int, T, T, T}}(undef, 0),  0, 0, false)
  end


  function AndersonAccelerator{T, BT, MT, RE}(dim::Int; mem::Int = 10, min_mem::Int = 3, λ::T = T(1e-8), activated::Bool = true) where {T <: AbstractFloat, RE <: AbstractRegularizer, BT <: AbstractBroydenType, MT <: AbstractMemory}
    mem, min_mem = arguments_check(mem, min_mem, dim, λ)

    F = zeros(T, dim, mem)
    X = zeros(T, dim, mem)
    G = zeros(T, dim, mem)
    M = zeros(T, mem, mem)
    Q = zeros(T, 0, 0)
    R = zeros(T, 0, 0)
    x_last = zeros(T, dim) 
    new(true, mem, min_mem, dim, 0, 0, zeros(Int,0), zeros(Int,0), zeros(Int,0), x_last, zeros(T, dim), zeros(T, dim),  zeros(T, dim), zeros(T, mem), F, X, G, M, Q, R, λ, activated, false, 0., 0., Vector{Tuple{Int, Symbol}}(undef, 0), Vector{Tuple{Int, T, T, T}}(undef, 0), 0, 0, false)
  end
 
  # specific constructor for Type 2 with QR decomposition
  function AndersonAccelerator{T, Type2{QRDecomp}, RestartedMemory, RE}(dim::Int; mem::Int = 10, min_mem::Int = 3, λ::T = T(1e-8), activated::Bool = true) where {T <: AbstractFloat, RE <: AbstractRegularizer}
    
    mem, min_mem = arguments_check(mem, min_mem, dim, λ)
    
    # for QR decomposition we don't need some of the caches
    x_last = zeros(T, 0)
    F = zeros(T, 0, 0)
    X = zeros(T, 0, 0)
    G = zeros(T, dim, mem)
    M = zeros(T, 0, 0)
    Q = zeros(T, dim, mem)
    R = zeros(T, mem, mem)
  
    new(true, mem, min_mem, dim, 0, 0, zeros(Int,0), zeros(Int,0), zeros(Int,0), x_last, zeros(T,dim), zeros(T, dim), zeros(T, dim),  zeros(T, dim), F, X, G, M, Q, R, λ, activated, false, 0., 0., Vector{Tuple{Int, Symbol}}(undef, 0), Vector{Tuple{Int, T, T, T}}(undef, 0), 0, 0, false)
  end

end

# define some default constructors for parameters
AndersonAccelerator(args...; kwargs...)  = AndersonAccelerator{Float64, Type2{QRDecomp}, RestartedMemory, NoRegularizer}(args...; kwargs...)
AndersonAccelerator{T}(args...; kwargs...) where {T <: AbstractFloat} = AndersonAccelerator{T,  Type2{QRDecomp}, RestartedMemory, NoRegularizer}(args...; kwargs...)
AndersonAccelerator{BT}(args...; kwargs...) where {BT <: AbstractBroydenType} = AndersonAccelerator{Float64, BT, RestartedMemory, NoRegularizer}(args...; kwargs...)
AndersonAccelerator{M}(args...; kwargs...) where {M <: AbstractMemory} = AndersonAccelerator{Float64, Type2{NormalEquations}, M, NoRegularizer}(args...; kwargs...)
AndersonAccelerator{RE}(args...; kwargs...) where {RE <: AbstractRegularizer} = AndersonAccelerator{Float64,  Type2{NormalEquations}, RestartedMemory, RE}(args...; kwargs...)

# restrict certain combination of parameters
AndersonAccelerator{T, Type2{QRDecomp}, RollingMemory, RE}(dim::Int; mem::Int = 10, min_mem::Int = 3, λ::T = T(1e-8), activated::Bool = true) where {T <: AbstractFloat, RE <: AbstractRegularizer}= throw(ArgumentError("Anderson - Type 2 with QR decomposition and RollingMemory not yet implemented."))
AndersonAccelerator{T, Type2{QRDecomp}, RestartedMemory, RE}(dim::Int; mem::Int = 10, min_mem::Int = 3, λ::T = T(1e-8), activated::Bool = true) where {T <: AbstractFloat, RE <: TikonovRegularizer} = throw(ArgumentError("Anderson - Type 2 with QR decomposition and regularisation not yet implemented."))

get_type(::AndersonAccelerator{T, BT, M, RE}) where {T, RE, BT, M} = BT
get_memory(::AndersonAccelerator{T, BT, M, RE}) where {T, RE, BT, M} = M

get_regularizer(::AndersonAccelerator{T, BT, M, RE}) where {T, RE, BT, M} = RE
is_active(aa::AndersonAccelerator) = aa.activated
was_successful(aa::AndersonAccelerator) = aa.success


function arguments_check(mem::Int, min_mem::Int, dim::Int, λ::T) where {T <: AbstractFloat}
  mem <= 1 && throw(ArgumentError("Memory has to be bigger than one."))
  min_mem <= 1 && throw(ArgumentError("Minimum memory has to be bigger than one."))
  dim <= 0 && throw(ArgumentError("Dimension has to be a positive integer."))
  (λ <= 0 || λ >= 1) && throw(ArgumentError("Regularisation parameter λ has to be between 0 and 1."))
  
  # mem shouldn't be bigger than the dimension
  mem = min(mem, dim)
  min_mem = min(min_mem, mem)

  return mem, min_mem
end

function log!(aa::AndersonAccelerator, iter::Int, status::Symbol)
  aa.activate_logging && push!(aa.acceleration_status, (num_iter, status))
  if status == :acc_guarded_declined
    aa.num_safe_declined += 1
  elseif status == :acc_guarded_accepted
    aa.num_accelerated_steps += 1
    aa.num_safe_accepted += 1
  elseif status == :acc_unguarded
    aa.num_accelerated_steps += 1   
  end
  return nothing
end


function activate!(aa::AndersonAccelerator) 
  aa.activated = true
end

function deactivate!(aa::AndersonAccelerator) 
  aa.activated = false
end


function restart!(aa::AndersonAccelerator{T, BT, MT, RE}) where {T <: AbstractFloat, RE <: AbstractRegularizer, BT <: AbstractBroydenType, MT <: AbstractMemory}
  empty_history!(aa)
end



function empty_history!(aa::AndersonAccelerator{T}) where {T <: AbstractFloat}
  # To save the deletion operation, we just reset the aa.iter pointer which 
  # tracks the number of filled columns 
  # aa.F .= 0;
  # aa.X .= 0;
  # aa.G .= 0;
  # aa.Q .= 0;
  # aa.R .= 0;

  # aa.f .= 0; we keep it for safeguarding
  aa.f_last .= 0;
  aa.g_last .= 0;
  aa.x_last .= 0;
  aa.eta .= 0;

  aa.iter = 0
  aa.init_phase = true
  aa.success = false
 
end

"Reset the pointer that tracks the number of valid cols in the cached history of past vectors."
function empty_caches!(aa::AndersonAccelerator)
  #to save time we can leave the old data in here, this is a potential source of problems if the indexing of valid data gets messed up
  # aa.F .= 0; 
  # aa.X .= 0;
  # aa.G .= 0;
  aa.iter = 0 #this is important as for the RestartedMemory holds information on how many rows are full
end

"""
  update!(aa, g, x)

- Update history of accelerator `aa` with iterates g = g(xi)
- Computes residuals f = x - g
"""
function update!(aa::AndersonAccelerator{T, BT, MT, RE}, g::AbstractVector{T}, x::AbstractVector{T}, num_iter::Int) where {T <: AbstractFloat, RE <: AbstractRegularizer, BT <: AbstractBroydenType, MT <: AbstractMemory}
  update_time_start = time()  
  if aa.activated
    # compute residual
    @. aa.f = x - g

    if aa.init_phase
      set_prev_vectors!(aa.x_last, aa.g_last, aa.f_last, x, g, aa.f)
      aa.init_phase = false
      return nothing
    end

    j = (aa.iter % aa.mem) + 1 # (aa.iter % aa.mem) number of cols filled, j is the next col where data should be entered

    if j == 1 && aa.iter != 0
      apply_memory_approach!(aa, num_iter) # for a RestartedMemory approach we want to flush the data cache matrices and start from scratch
    end

    # store Δx, Δg, Δf in X, G, F
    fill_caches!(j, aa.X, aa.G, aa.F, x, g, aa.f, aa.x_last, aa.g_last, aa.f_last)

    # set previous values for next iteration
    set_prev_vectors!(aa.x_last, aa.g_last, aa.f_last, x, g, aa.f)

    aa.iter += 1
    aa.update_time += time() - update_time_start
  end
  return nothing
end

# This method is a copy of the method above to test TypeII with QR decomposition
function update!(aa::AndersonAccelerator{T, Type2{QRDecomp}, RestartedMemory, RE}, g::AbstractVector{T}, x::AbstractVector{T}, num_iter::Int) where {T <: AbstractFloat, RE <: AbstractRegularizer}
  if aa.activated
    update_time_start = time()  

    # compute residual
    @. aa.f = x - g

    if aa.init_phase
      set_prev_vectors!(aa.g_last, aa.f_last, g, aa.f)
      aa.init_phase = false
      return nothing
    end

    j = (aa.iter % aa.mem) + 1 # (aa.iter % aa.mem) number of cols filled, j is the next col where data should be entered

    if j == 1 && aa.iter != 0
      apply_memory_approach!(aa, num_iter) # for a RestartedMemory approach we want to flush the data cache matrices and start from scratch
    end

    #G[:, j] = Δg
    fill_delta!(j, aa.G, g, aa.g_last)
    
    # use f_last memory to store Δf = f - f_last 
    compute_Δf!(aa.f_last, aa.f)

    # QR decomposition step
    qr!(aa.Q, aa.R, aa.f_last, j)

    # set previous values for next iteration
    set_prev_vectors!(aa.g_last, aa.f_last, g, aa.f)

    aa.iter += 1
    aa.update_time += time() - update_time_start
  end
  return nothing
end

# adapted from https://github.com/JuliaNLSolvers/NLsolve.jl
# add new column to QR-Factorisation: F = QR 
function qr!(Q::AbstractMatrix{T}, R::AbstractMatrix{T}, Δf::AbstractVector{T}, j::Int) where {T <: AbstractFloat}
  
  n, m = size(Q)
  n == length(Δf) || throw(DimensionMismatch())
  m == LinearAlgebra.checksquare(R) || throw(DimensionMismatch())
  1 ≤ j ≤ m || throw(ArgumentError())

  @inbounds for k in 1:(j - 1)
    qk = uview(Q, :, k)
    ip = dot(qk, Δf)
    R[k, j] = ip

    axpy!(-ip, qk, Δf)
  end
  @inbounds begin
    nrm_f = norm(Δf, 2)
    R[j, j] = nrm_f
    @. Q[:, j] = Δf / nrm_f
  end
  return nothing
end


" Update the history matrices X = [Δxi, Δxi+1, ...], G = [Δgi, Δgi+1, ...] and F = [Δfi, Δfi+1, ...]."
function fill_caches!(j::Int, X::AbstractMatrix{T}, G::AbstractMatrix{T}, F::AbstractMatrix{T}, x::AbstractVector{T}, g::AbstractVector{T}, f::AbstractVector{T}, x_last::AbstractVector{T}, g_last::AbstractVector{T}, f_last::AbstractVector{T}) where {T <: AbstractFloat}
  # fill memory with deltas
  @. X[:, j] = x - x_last # Δx
  @. G[:, j] = g - g_last # Δg
  @. F[:, j] = f - f_last # Δf
  return nothing
end

" Update a single history matrices V = [vgi, vgi+1, ...] ."
function fill_delta!(j::Int, V::AbstractMatrix{T}, v::AbstractVector{T}, v_last::AbstractVector{T}) where {T <: AbstractFloat}
  @inbounds for i in eachindex(v)
    V[i, j] = v[i] - v_last[i]
  end
end

"Compute Δf = f - f_last and store result in f_last."
function compute_Δf!(f_last::AbstractVector{T}, f::AbstractVector{T}) where {T <: AbstractFloat}
  
  @inbounds for i in eachindex(f_last)
    f_last[i] *= -one(T)
    f_last[i] += f[i]
  end
end

"Store a copy of the last x, g, f to be able to compute Δx, Δg, Δf at the next step."
function set_prev_vectors!(x_last::AbstractVector{T}, g_last::AbstractVector{T}, f_last::AbstractVector{T}, x::AbstractVector{T}, g::AbstractVector{T}, f::AbstractVector{T}) where {T <: AbstractFloat}
  @. x_last = x
  @. g_last = g
  @. f_last = f
  return nothing
end

function set_prev_vectors!(g_last::AbstractVector{T}, f_last::AbstractVector{T}, g::AbstractVector{T}, f::AbstractVector{T}) where {T <: AbstractFloat}
  @. g_last = g
  @. f_last = f
  return nothing
end

"Depending on the AndersonAccelerator parameter AbstractMemory, dispatch on the correct method that handles the case when the memory buffers are full."
function apply_memory_approach!(aa::AndersonAccelerator{T, BT, RestartedMemory, RE}, num_iter) where {T, RE, BT}
    empty_caches!(aa)
    aa.activate_logging && log!(aa, num_iter, :memory_full)
    return true
end
apply_memory_approach!(aa::AndersonAccelerator{T, BT, RollingMemory, RE}, num_iter) where {T, RE, BT} = nothing

"Anderson type-1: eta = M^{-1} X' where M = (X'F)."
function assemble_inv_matrix!(M::AbstractMatrix{T}, X::AbstractMatrix{T}, F::AbstractMatrix{T}, aa::AndersonAccelerator{T, Type1, ME, RE}) where {T <: AbstractFloat, RE <: AbstractRegularizer, ME <: AbstractMemory}
  mul!(M, X', F)
end
"Anderson type-2: eta = M^{-1} X' where M = (F'F)."
function assemble_inv_matrix!(M::AbstractMatrix{T}, X::AbstractMatrix{T}, F::AbstractMatrix{T}, aa::AndersonAccelerator{T, Type2{NormalEquations}, ME, RE}) where {T <: AbstractFloat, RE <: AbstractRegularizer, ME <: AbstractMemory}
  mul!(M, F', F)
end

"Recombine past iterates to compute an accelerated point. Overwrite g with accelerated point."
function accelerate!(g::AbstractVector{T}, x::AbstractVector{T}, aa::AndersonAccelerator{T, BT, ME, RE}, num_iter::Int) where {T <: AbstractFloat, BT <: AbstractBroydenType, ME <: AbstractMemory, RE <: AbstractRegularizer}
  aa.success = false
  l = min(aa.iter, aa.mem) #number of columns filled with data
  if l < aa.min_mem
     aa.activate_logging && push!(aa.acceleration_status, (num_iter, :not_enough_cols))
     return nothing
  end
  accelerate_time_start = time()

  eta = uview(aa.eta, 1:l)
  X = uview(aa.X, :, 1:l)
  M = uview(aa.M, 1:l, 1:l)
  G = uview(aa.G, :, 1:l)
  F = uview(aa.F, :, 1:l)

  # compute eta depending on method type
  # Type-1: eta = (X'F)^{-1}X'f
  # Type-2: eta = (F'F)^{-1}F'f
  assemble_inv_matrix!(M, X, F, aa)
  initialise_eta!(eta, aa, X, F)
  # aa.eta = aa.M  \ eta 
  info = solve_linear_sys!(M, eta, aa)

  if (info < 0 || norm(eta, 2) > 1e4)
    if info < 0
      aa.activate_logging && push!(aa.acceleration_status, (num_iter, :fail_singular))
    elseif norm(eta, 2) > 1e4
      aa.activate_logging && push!(aa.acceleration_status, (num_iter, :fail_eta_norm))
    end
    aa.accelerate_time += time() - accelerate_time_start
    return nothing
  else
    # calculate the accelerated candidate point and overwrite g
    # g = g - G * eta
    BLAS.gemv!('N', -one(T), G, eta, one(T), g)
    aa.accelerate_time += time() - accelerate_time_start
    aa.success = true
    return nothing
  end
end

"Recombine past iterates to compute an accelerated point. Overwrite g with accelerated point. Uses QR-decomposition."
function accelerate!(g::AbstractVector{T}, x::AbstractVector{T}, aa::AndersonAccelerator{T, Type2{QRDecomp}, RestartedMemory, RE}, num_iter::Int) where {T <: AbstractFloat, RE <: AbstractRegularizer}
  aa.success = false
  l = min(aa.iter, aa.mem) #number of columns filled with data
  if l < aa.min_mem
     aa.activate_logging && push!(aa.acceleration_status, (num_iter, :not_enough_cols))
     return nothing
  end
  accelerate_time_start = time()

  eta = uview(aa.eta, 1:l)
  G = uview(aa.G, :, 1:l)
  Q = uview(aa.Q, :, 1:l)
  R = uview(aa.R, 1:l, 1:l) 

  # solve least squares problem ||f_k - η Fk ||_2 where Fk = QR 
  # initialise_eta!(eta, aa, X, F)
  mul!(eta, Q', aa.f) 
  info = solve_linear_sys!(R, eta, aa)
   if info < 0
    aa.activate_logging && push!(aa.acceleration_status, (num_iter, :fail_singular))
    aa.accelerate_time += time() - accelerate_time_start
    return nothing
  end
  
  # TODO: maybe replace this with a check of the condition number of R
  if norm(eta, 2) > 1e4
    aa.activate_logging && push!(aa.acceleration_status, (num_iter, :fail_eta_norm))
    aa.accelerate_time += time() - accelerate_time_start
    return nothing
  else
    # calculate the accelerated candidate point
    # g = g - G * eta
    BLAS.gemv!('N', -one(T), G, eta, one(T), g)
        
    aa.accelerate_time += time() - accelerate_time_start
    aa.success = true
    return nothing
  end
end


"Anderson Type 1: Initialise eta = X'f."
function initialise_eta!(eta::AbstractVector{T}, aa::AndersonAccelerator{T, Type1}, X::AbstractMatrix{T}, F::AbstractMatrix{T}) where {T <: AbstractFloat}
  mul!(eta, X', aa.f)
end

"Anderson Type 2: Initialise eta = F'f."
function initialise_eta!(eta::AbstractVector{T}, aa::AndersonAccelerator{T, Type2{NormalEquations}}, X::AbstractMatrix{T}, F::AbstractMatrix{T}) where {T <: AbstractFloat}
  mul!(eta, F', aa.f)
end



# BLAS gesv! wrapper with error handling
# solve A X = B and for X and save result in B
function _gesv!(A, B)
  try
    LinearAlgebra.LAPACK.gesv!(A, B)
    return 1
    catch
      return -1
  end
end

  
function solve_linear_sys!(M::AbstractMatrix{T}, eta::AbstractVector{T}, aa::AndersonAccelerator{T, BT, ME, TikonovRegularizer}) where {T <: AbstractFloat, BT <: AbstractBroydenType, ME <: AbstractMemory}
   # add regularisation term
  for i = 1:size(M, 1)
    M[i, i] += β
  end
  # solve regularised problem
  info = _gesv!(M, eta)
end

function solve_linear_sys!(M::AbstractMatrix{T}, eta::AbstractVector{T}, aa::AndersonAccelerator{T, BT, ME, NoRegularizer}) where {T <: AbstractFloat, BT <: AbstractBroydenType, ME <: AbstractMemory}
  info = _gesv!(M, eta)
end


function solve_linear_sys!(R::AbstractMatrix{T}, eta::AbstractVector{T}, aa::AndersonAccelerator{T, Type2{QRDecomp}, ME, NoRegularizer}) where {T <: AbstractFloat, ME <: AbstractMemory}
  try
    LinearAlgebra.ldiv!(R, eta)
    # LinearAlgebra.BLAS.trsv!('U', 'N', 'N', R, eta) # seems to be equally fast, but can lead to NaN in eta without errors in openblas
    return 1
   catch
     return -1
  end
end
