module COSMOAccelerators

using SparseArrays, LinearAlgebra
export update!, accelerate!, restart!, log!, was_successful 
    include("./abstract_accelerator.jl")
    include("./empty_accelerator.jl")
    include("./anderson_accelerator.jl")


    function version()
        v"0.1.0"
    end
end #end module
