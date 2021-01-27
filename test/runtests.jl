using COSMOAccelerators
using Test

@testset "COSMOAccelerators.jl" begin
    include("proximal_operators.jl")
    include("empty.jl")
    include("anderson.jl")
end
