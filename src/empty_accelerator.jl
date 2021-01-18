export EmptyAccelerator
# ---------------------------
# EmptyAccelerator
# ---------------------------

"""
    EmptyAccelerator{T} <: AbstractAccelerator

Accelerator that does nothing.
"""
struct EmptyAccelerator{T} <: AbstractAccelerator 
  function EmptyAccelerator{T}() where {T <: AbstractFloat}
      return new{T}()
  end
end
EmptyAccelerator(args...; kwargs...) = EmptyAccelerator{Float64}(args...; kwargs...)
EmptyAccelerator{T}(dim::Int64) where {T <: AbstractFloat} = EmptyAccelerator{T}()

function update!(ea::EmptyAccelerator{T}, args...; kwargs...) where {T <: AbstractFloat}
  return nothing
end

function accelerate!(g::AbstractVector{T}, x::AbstractVector{T}, ea::EmptyAccelerator{T}, args...; kwargs... ) where {T <: AbstractFloat}
  return false
end
restart!(ea::EmptyAccelerator, args...; kwargs...) = nothing
activate!(ea::EmptyAccelerator, args...; kwargs...) = nothing
log!(ea, args...; kwargs...) = nothing
was_successful(ea::EmptyAccelerator) = false
is_active(ea::EmptyAccelerator) = false
get_memory_size(ea::EmptyAccelerator) = 0