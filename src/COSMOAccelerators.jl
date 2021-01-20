module COSMOAccelerators

using SparseArrays, LinearAlgebra, UnsafeArrays
export update!, accelerate!, restart!, activate!, log!, was_successful, is_active
    include("./abstract_accelerator.jl")
    include("./empty_accelerator.jl")
    include("./anderson_accelerator.jl")


    function version()
        v"0.1.0"
    end
end #end module
