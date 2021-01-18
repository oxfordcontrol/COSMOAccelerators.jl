using COSMOAccelerators

@testset "EmptyAccelerator" begin
    ea = EmptyAccelerator()
    @test typeof(ea) == EmptyAccelerator{Float64}
    iter = 1
    @test update!(ea, zeros(2), zeros(2),iter) == nothing
    @test accelerate!(zeros(2), zeros(2), ea, iter) == false
    
    @test get_memory_size(ea) == 0
    @test was_successful(ea) == false
    @test restart!(ea) == nothing
    @test is_active(ea) == false

end