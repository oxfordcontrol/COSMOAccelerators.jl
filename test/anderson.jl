using COSMOAccelerators, Test, Random, LinearAlgebra, SparseArrays

@testset "AndersonAccelerator" begin

    rng = Random.MersenneTwister(2)



    @testset "Constructors" begin
        dim = 10
        # default: Float64, Type 2, QR decomposition, NoRegularizer, RestartedMemory
        aa = AndersonAccelerator(dim, mem = 10)
        @test typeof(aa) == AndersonAccelerator{Float64, Type2{QRDecomp}, RestartedMemory, NoRegularizer}
        
        @test get_memory_size(aa) == 10
        @test is_active(aa) == true

        aa = AndersonAccelerator{Float32}(dim)
        @test typeof(aa) == AndersonAccelerator{Float32, Type2{QRDecomp}, RestartedMemory, NoRegularizer} 

        aa = AndersonAccelerator{Type1}(dim)
        @test typeof(aa) == AndersonAccelerator{Float64, Type1, RestartedMemory, NoRegularizer} 

        aa = AndersonAccelerator{TikonovRegularizer}(dim)
        @test typeof(aa) == AndersonAccelerator{Float64,  Type2{NormalEquations}, RestartedMemory, TikonovRegularizer} 

        aa = AndersonAccelerator{RollingMemory}(dim)
        @test typeof(aa) == AndersonAccelerator{Float64,  Type2{NormalEquations}, RollingMemory, NoRegularizer} 

        @test_throws ArgumentError AndersonAccelerator(dim, mem = 0)
        @test_throws ArgumentError AndersonAccelerator(dim, min_mem = 0)
        @test_throws ArgumentError AndersonAccelerator(dim, λ = -1.)
        @test_throws ArgumentError AndersonAccelerator{Float64, Type2{QRDecomp}, RollingMemory, NoRegularizer}(dim)
        @test_throws ArgumentError AndersonAccelerator{Float64, Type2{QRDecomp}, RestartedMemory, TikonovRegularizer}(dim)
    end




    @testset "Anderson Acceleration (RollingMemory)" begin
        # add a history of vectors
        n = 2
        m = 4
        hist_length = 10
        G = randn(rng, n + m, hist_length)
        X = randn(rng, n + m, hist_length)
        F = X - G
        mem = 5
        aa = AndersonAccelerator{Float64, Type1, RollingMemory, NoRegularizer}(m + n, mem = mem);

        # add more vectors than memory, it should start appending at the beginning of the matrix caches again
        for k = 1:mem+2
            update_history!(aa, G[:, k], X[:, k], k) 
        end 
        X_ref = X[:,7] - X[:, 6]
        F_ref = F[:,7] - F[:, 6]
        G_ref = G[:,7] - G[:, 6]
        @test aa.X[:, 1] == X_ref
        @test aa.G[:, 1] == G_ref
        @test aa.F[:, 1] == F_ref
        return nothing

    end


    @testset "Anderson Acceleration (RestartedMemory)" begin

        # add a history of vectors
        n = 2
        m = 4
        hist_length = 10
        G = randn(rng, n + m, hist_length)
        X = randn(rng, n + m, hist_length)
        F = X - G
        mem = 5
        aa = AndersonAccelerator{Float64, Type2{NormalEquations}, RestartedMemory, NoRegularizer}(m + n, mem = mem);

        # lets add a bunch of inputs and double check that everything is stored correctly
        @test aa.init_phase

        # add first vector pair
        update_history!(aa, G[:, 1], X[:, 1], 1)
        @test !aa.init_phase 
        @test aa.iter == 0
        @test aa.x_last == X[:, 1]
        @test aa.g_last == G[:, 1]
        @test aa.f_last == X[:, 1] - G[:, 1]
        @test iszero(aa.F) && iszero(aa.G) && iszero(aa.F) && iszero(aa.M)
        
        # add second vector
        update_history!(aa, G[:, 2], X[:, 2], 2) 
        @test aa.iter == 1
        @test aa.F[:, 1] ==  F[:, 2] - F[:, 1] 
        @test aa.X[:, 1] ==  X[:, 2] - X[:, 1] 
        @test aa.G[:, 1] ==  G[:, 2] - G[:, 1] 

        # add third and fourth vector pairs
        k = 3
        update_history!(aa, G[:, k], X[:, k], 3) 
        k += 1
        update_history!(aa, G[:, k], X[:, k], 4) 
        @test aa.iter == 3
        # let's check the matrices
        X_ref = [X[:,2] - X[:, 1] X[:, 3] - X[:, 2] X[:, 4] - X[:, 3] ]   
        F_ref = [F[:,2] - F[:, 1] F[:, 3] - F[:, 2] F[:, 4] - F[:, 3] ]   
        G_ref = [G[:,2] - G[:, 1] G[:, 3] - G[:, 2] G[:, 4] - G[:, 3] ]   
        @test aa.X[:,1:3] ==  X_ref
        @test aa.F[:,1:3] ==  F_ref
        @test aa.G[:,1:3] ==  G_ref

        # let's test the accelerated vector
        v = copy(G[:, k])
        v_prev = copy(X[:, k])
        f = F[:, k]
        eta = F_ref' * f  
        eta = F_ref' * F_ref \ eta 
        v_acc_ref = v - (X_ref - F_ref) * eta 
        accelerate!(v, v_prev, aa, k)
        @test isapprox(v, v_acc_ref, atol = 1e-8)
        

        # let's add more vectors to exceed memory and check that memory gets flushed 
        k += 1
        update_history!(aa, G[:, k], X[:, k], k) 
        k += 1 
        update_history!(aa, G[:, k], X[:, k], k) 
        # memory should now be full
        @test aa.iter % aa.mem == 0

        # add one more and double check it is added in first position
        k += 1
        update_history!(aa, G[:, k], X[:, k], k) 
 
        
        @test aa.X[:, 1] == X[:,7] - X[:, 6] 
        @test aa.F[:, 1] == F[:,7] - F[:, 6] 
        @test aa.G[:, 1] == G[:,7] - G[:, 6] 
        @test aa.iter == 1
        @test aa.init_phase == false
        return nothing


    end

    @testset "Gram-Schmidt Orthogonalisation method" begin
        F = LinearAlgebra.qr(rand(rng, 8, 4))    
        V = Matrix(F.Q)

        # now test our method
        V2 = [V zeros(8)]

        COSMOAccelerators.qr!(V2, zeros(5, 5), v2, 5)
        v2 = V2[:, 5]
        # is new vector orthogonal to all other vectors in V
        for k = 1:4
            @test dot(V[:, k], v2) <= 1e-12
        end 
        @test isapprox(norm(v2), 1., atol = 1e-8)
        return nothing


    end

end