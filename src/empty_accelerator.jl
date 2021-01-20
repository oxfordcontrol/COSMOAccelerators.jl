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


function update!(ea::EmptyAccelerator, args...; kwargs...)
  return nothing
end

function accelerate!(g::AbstractVector{T}, x::AbstractVector{T}, ea::EmptyAccelerator, args...; kwargs... ) where {T <: AbstractFloat}
  return false
end
restart!(ea::EmptyAccelerator, args...; kwargs...) = nothing
log!(ea, args...; kwargs...) = nothing
was_successful(ea::EmptyAccelerator) = false
