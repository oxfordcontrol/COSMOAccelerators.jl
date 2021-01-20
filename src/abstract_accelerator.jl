# An abstract type for fixed-point acceleration methods
# Fixed point problem x = g(x) with residual f(x) = x - g(x)
export AbstractAccelerator
"""
    AbstractAccelerator

Abstract supertype for acceleration objects that can be used to speed up a fixed-point iterations g = g(x) of a nonexpansive operator `g`. 
They must implement the following methods to communicate with the fixed-point algorithm:
  - update!(aa::AbstractAccelerator, g, x, num_iter) #stores the fixed-point iterates
  - accelerate!(g::AbstractVector, x, aa::AbstractAccelerator, num_iter) #recombines past iterates to determine an accelerated point and overwrites `g`
  - restart!(aa::AbstractAccelerator, args...; kwargs...) # algorithm tells the accelerator to restart

The algorithm has to be able to query the following information:
  - was_successful(aa::AbstractAccelerator) --> Bool #indicate whether accelerate! was succesful at the last iteration

The following method is optional:
  - log!(aa::AbstractAccelerator, args...; kwargs...) # algorithm tells accelerator to log certain information for debugging
"""
abstract type AbstractAccelerator end

# update!(aa::AbstractAccelerator, g, x, iter) = nothing
# accelerate!(g, x, aa::AbstractAccelerator, iter) = nothing
# was_successful(aa::AbstractAccelerator) = false
# restart!(aa::AbstractAccelerator) = nothing

log!(aa::AbstractAccelerator, args...; kwargs...) = false
