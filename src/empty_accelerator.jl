export EmptyAccelerator
# ---------------------------
# EmptyAccelerator
# ---------------------------

"""
    EmptyAccelerator <: AbstractAccelerator

An accelerator that does nothing.
"""
struct EmptyAccelerator <: AbstractAccelerator 
  function EmptyAccelerator() 
      return new()
  end
end
EmptyAccelerator(args...; kwargs...) = EmptyAccelerator()


function update!(ea::EmptyAccelerator, args...)
  return nothing
end

function accelerate!(g::AbstractVector{T}, x::AbstractVector{T}, ea::EmptyAccelerator, args... ) where {T <: AbstractFloat}
  return false
end
restart!(ea::EmptyAccelerator) = nothing
log!(ea, args...; kwargs...) = nothing
was_successful(ea::EmptyAccelerator) = false
