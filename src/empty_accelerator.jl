export EmptyAccelerator
# ---------------------------
# EmptyAccelerator
# ---------------------------

"""
    EmptyAccelerator <: AbstractAccelerator

Accelerator that does nothing.
"""
struct EmptyAccelerator <: AbstractAccelerator 
  function EmptyAccelerator() 
      return new()
  end
end
EmptyAccelerator(args...; kwargs...) = EmptyAccelerator()


function update!(ea::EmptyAccelerator, args...; kwargs...)
  return nothing
end

function accelerate!(g::AbstractVector{T}, x::AbstractVector{T}, ea::EmptyAccelerator, args...; kwargs... ) where {T <: AbstractFloat}
  return false
end
restart!(ea::EmptyAccelerator, args...; kwargs...) = nothing
activate!(ea::EmptyAccelerator, args...; kwargs...) = nothing
log!(ea, args...; kwargs...) = nothing
was_successful(ea::EmptyAccelerator) = false
is_active(ea::EmptyAccelerator) = false
get_memory_size(ea::EmptyAccelerator) = 0