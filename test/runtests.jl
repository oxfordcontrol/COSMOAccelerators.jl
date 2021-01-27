using COSMOAccelerators
using Test

@testset "COSMOAccelerators.jl" begin
    include("proximal_operators.jl")
    include("empty.jl")
    include("anderson.jl")
end

using UnsafeArrays
R = rand(5, 5)

eta = rand(10)

Rv = uview(R, 1:4, 1:4)
etav = uview(eta, 1:4)

@code_llvm LinearAlgebra.ldiv!(Rv, etav


@which LAPACK.trtrs!('U', 'N', 'N', Rv, etav)