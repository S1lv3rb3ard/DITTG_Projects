using FFTW
using TimerOutputs
using Distributed
using LinearAlgebra
@everywhere using SharedArrays

@everywhere struct Trilayer
    # Unrotated monolayer parameters
    E::Array{Float64,2}
    P::Array{Float64,1}

    # Twist, rotation
    θ::Array{Float64,1}
    R::Array{Array{Float64,2},1}
    iR::Array{Array{Float64,2},1}

    # Twisted layers parameterization
    tE::Array{Array{Float64,2},1}
    invtE::Array{Array{Float64,2},1}
    tP::Array{Array{Float64,1},1}

    # Effective elastic tensor in hull space
    C::Array{Float64,2}

    function Trilayer(E, P, θ1, θ3, K, G)

        ## Symmetries and rotations:
        R = [   [cos(θ1) sin(θ1); -sin(θ1) cos(θ1)],
                [1       0;         0         1   ],
                [cos(θ3) sin(θ3); -sin(θ3) cos(θ3)] ];

        ## Lattice parameters for twisted layers:
        θ = [θ1,0.,θ3]
        tE = R .* [E]
        tP = R .* [P]

        ## Intra-layer elastic tensor, with stress tensor linear ordering [E_11; E_12; E_21; E_22]
        C = [K+G 0 0 K-G ; 0 G G 0; 0 G G 0; K-G 0 0 K+G]

        new(E, P, θ, R, inv.(R), tE, 2*pi*inv.(tE), tP, C)
    end
end

# Discretized Hull structure (same as the bilayer case)
struct Hull
    tl::Trilayer
    N::Int64
    hN::Int64
    Γ0::SharedArray{Float64,2}
    tK::Array{ComplexF64,6}
    Hess_Elastic::SharedArray{ComplexF64,3}
    Precon_Elastic::Array{ComplexF64,3}
    plan
    iplan
    permutations::SharedArray{Int64,2}

    cache_u::SharedArray{ComplexF64,6}
    cache_grad::SharedArray{ComplexF64,6}
    cache_ur::SharedArray{Float64,6}
    cache_misfit::SharedArray{Float64,1}
    cache_grad_misfit::SharedArray{Float64,3}
    cache_shifts::SharedArray{Float64,2}
    cache_gsfe_computation::SharedArray{Float64,2}

    cache_Misfits::SharedArray{Float64,1}

    # Generate Hull discretization
    function Hull(tl::Trilayer, N)
        s =        (0:N-1)            .* ones(Int,1,N,N,N);
        t = reshape(0:N-1, (1,N,1,1)) .* ones(Int,N,1,N,N);
        v = reshape(0:N-1, (1,1,N,1)) .* ones(Int,N,N,1,N);
        w = reshape(0:N-1, (1,1,1,N)) .* ones(Int,N,N,N,1);

        Γ0    = (2*pi/N) .* (hcat(s[:], t[:])')

        # Prepare the arrays encoding the permutations, u2^c(γ_2,γ3,0) and u2^c(γ1,γ2,0)
        arr = hcat(-s[:], -t[:], v[:]-s[:], w[:]-t[:])'
        arr[arr .< 0] .+= N

        # mapping values to index
        ind1 = (arr[1, :] .+ 1) .+ arr[2, :] * N .+ arr[3, :] * N^2 .+ arr[4, :] * N^3
        ind2 = (arr[3, :] .+ 1) .+ arr[4, :] * N .+ arr[1, :] * N^2 .+ arr[2, :] * N^3
        permutations = cat(ind1[:], ind2[:]; dims = 2)

        # Matrices encoding the real space gradient directions in configuration space
        T1 = cat(   (inv(tl.tE[1]) - inv(tl.tE[2]))',
                    (inv(tl.tE[2]) - inv(tl.tE[1]))',
                    (inv(tl.tE[3]) - inv(tl.tE[2]))',   dims = 7)
        T2 = cat(   (inv(tl.tE[1]) - inv(tl.tE[3]))',
                    (inv(tl.tE[2]) - inv(tl.tE[3]))',
                    (inv(tl.tE[3]) - inv(tl.tE[1]))',   dims = 7)

        # Fourier multipliers for the real space gradient

        hN = div(N, 2) + 1
        k = 2*pi*[0:hN-1;hN-N:-1]
        K = cat(   reshape(k[1:hN], (1,1,hN,1,1,1,1)) .* ones(1,1,1,N,1,1,1),
                   reshape(k, (1,1,1,N,1,1,1)) .* ones(1,1,hN,1,1,1,1);  dims = 2)
        tK1 = sum(T1 .* K; dims = 2)
        K = cat(   reshape(k, (1,1,1,1,N,1,1)) .* ones(1,1,1,1,1,N,1),
                   reshape(k, (1,1,1,1,1,N,1)) .* ones(1,1,1,1,N,1,1);  dims = 2)
        tK2 = sum(T2 .* K; dims = 2)

        tK = reshape(tK1 .+ tK2, (2,hN,N,N,N,3))

        # Pointwise outer product
        tKtK =  reshape(tK, (2,1,hN,N,N,N,3)) .* reshape(tK, (1,2,hN,N,N,N,3))
        # Hessian matrix (diagonal in Fourier space with i,j,k,l indices, layer index p)
        C = reshape(tl.C,(2,2,2,2))
        Hess = zeros(ComplexF64, 2,2,hN*N^3*3)
        for d=1:2, b=1:2
            Hess .+=    reshape(C[:,b,:,d],          (2,2,1)) .*
                        reshape(tKtK[b,d,:,:,:,:,:], (1,1,hN*N^3*3))
        end
        Precon = zeros(ComplexF64, 2,2,hN*N^3*3)
        for j=1:3
            for i=2:hN*N^3
                idx = i+(j-1)*hN*N^3
                T = svd(Hess[:,:,idx])
                Precon[:,:,idx] = T.U * Diagonal( 1. ./ sqrt.(0.01 .+ T.S.^2) ) * T.Vt
            end
        end
        Precon = reshape(Precon, (2,2,hN*N^3*3))

        # allocate memory for FFTs
        A = zeros(ComplexF64, (2,hN,N,N,N,3))
        iplan = 1/N^2 * plan_brfft(A, N, (2,3,4,5), flags=FFTW.MEASURE)     # to apply ifft, use plan*P -> result is the ifft of A.
        A = zeros(Float64, (2,N,N,N,N,3))
        plan = 1/N^2 * plan_rfft(A, (2,3,4,5), flags=FFTW.MEASURE)          # to apply fft, use iplan*P -> result is the fft of A.

        cache_u = zeros(ComplexF64, (2,hN,N,N,N,3))
        cache_grad = zeros(ComplexF64, (2,hN,N,N,N,3))
        cache_ur = zeros(Float64, (2,N,N,N,N,3))
        cache_misfit = zeros(Float64, N^4)
        cache_grad_misfit = zeros(Float64, (2,N^4,3))
        cache_shifts = zeros(Float64, (2,N^4))
        cache_gsfe_computation = zeros(Float64, (2,N^4))

        cache_Elastics = SharedArray{Float64,1}(nworkers())
        cache_Misfits = SharedArray{Float64,1}(nworkers())

        new(tl, N, hN, Γ0, tK, Hess, Precon, plan, iplan, permutations, cache_u, cache_grad, cache_ur, cache_misfit, cache_grad_misfit, cache_shifts, cache_gsfe_computation, cache_Misfits)
    end
end

# Real-space gradient
function ∇(u::Union{Array{Float64, 2},Array{Float64, 6}}, hull::Hull)
    uf = hull.plan * reshape(u, (2, N, N, N, N, 3))
    ux = hull.iplan * (1im .* hull.tK[1:1,:,:,:,:,:] .* uf)
    uy = hull.iplan * (1im .* hull.tK[2:2,:,:,:,:,:] .* uf)
    return [[ux[1,:,:,:,:,:]] [uy[1,:,:,:,:,:]];
            [ux[2,:,:,:,:,:]] [uy[2,:,:,:,:,:]]]
end

function ∇(u::Array{ComplexF64, 6}, hull::Hull)
    ux = hull.iplan * (1im .* hull.tK[1:1,:,:,:,:,:] .* u)
    uy = hull.iplan * (1im .* hull.tK[2:2,:,:,:,:,:] .* u)
    return [[ux[1,:,:,:,:,:]] [uy[1,:,:,:,:,:]];
            [ux[2,:,:,:,:,:]] [uy[2,:,:,:,:,:]]]
end

# Preconditioning utilities
function LinearAlgebra.ldiv!(out::Array{ComplexF64,6}, P::Array{ComplexF64,3}, A::Array{ComplexF64,6})
    len = prod(size(out)[2:6])
    vOut = reshape(out, 2, len)
    vA = reshape(A, 2, len)
    @timeit to "ldiv!" @inbounds @fastmath @simd for i=1:len
        vOut[1,i] = P[1,1,i] * vA[1,i] + P[2,1,i] * vA[2,i]
        vOut[2,i] = P[1,2,i] * vA[1,i] + P[2,2,i] * vA[2,i]
    end
    # @sync begin
    #     len = size(P,3)
    #     for pw in workers()
    #         @async remotecall_wait(g_elastic_shared!, pw,
    #                     reshape(out, (2,len)), len,
    #                     P, reshape(A, (2,len)))
    #     end
    # end
end
function LinearAlgebra.dot(A::Array{ComplexF64,6}, B::Array{ComplexF64,6})
    inc = prod(size(A)[1:2])
    len = prod(size(A)[3:end])
    @timeit to "dot" return real(
                           2.0 *  BLAS.dotc(len*inc, pointer(A), 1, pointer(B), 1)
                                - BLAS.dotc(len, pointer(A,1), inc, pointer(B,1), inc)
                                - BLAS.dotc(len, pointer(A,2), inc, pointer(B,2), inc)
                                )
end


@everywhere function myrange(idx::Int64, M::Int64)
    nwks = length(workers())
    if (idx < nwks)
        return rge = 1 + div((idx-1)*M, nwks) : div(idx*M, nwks)
    else # idx = nwks
        return rge = 1 + div((idx-1)*M, nwks) : M
    end
end

@everywhere function g_elastic_shared!(  storage_g::SharedArray{ComplexF64,2},
                                         M::Int64,
                                         H::SharedArray{ComplexF64,3},
                                         U::SharedArray{ComplexF64,2})
    idx = indexpids(storage_g)
    rge = myrange(idx, M)

    @inbounds @fastmath @simd for i=rge
        storage_g[1,i] = H[1,1,i] * U[1,i] + H[2,1,i] * U[2,i]
        storage_g[2,i] = H[1,2,i] * U[1,i] + H[2,2,i] * U[2,i]
    end
    nothing
end

@everywhere function gsfe_kernel( rge,
                                  f::Union{Nothing, Float64},
                                  g::Union{Nothing, SubArray{Float64,2}},
                                  invEi::Array{Float64,2},
                                  Γ0::SharedArray{Float64,2},
                                  uj::SubArray{Float64,2},
                                  u2::SubArray{Float64,2},
                                  permutation::SubArray{Int64,1},
                                  misfit::SharedArray{Float64,1},
                                  shifts::SharedArray{Float64,2},
                                  g_misfit::SharedArray{Float64,2})

@inbounds @fastmath @simd for i=rge
        shifts[1,i] = Γ0[1,i] + invEi[1,1]*(u2[1, permutation[i]] - uj[1,i]) +
                                      invEi[1,2]*(u2[2, permutation[i]] - uj[2,i])
        shifts[2,i] = Γ0[2,i] + invEi[2,1]*(u2[1, permutation[i]] - uj[1,i]) +
                                      invEi[2,2]*(u2[2, permutation[i]] - uj[2,i])
    end
    if f != nothing
@views  GSFE!(misfit[rge], shifts[:,rge])
        @inbounds for i=rge
            f += misfit[i]
        end
    end
    if g != nothing
@views  gradient_GSFE!(g_misfit[:,rge], shifts[:,rge])
        @inbounds @fastmath @simd for i=rge
            g[1,i] -= 0.5 * (   invEi[1,1] * g_misfit[1,i]
                              + invEi[2,1] * g_misfit[2,i])
            g[2,i] -= 0.5 * (   invEi[1,2] * g_misfit[1,i]
                              + invEi[2,2] * g_misfit[2,i])
        end
    end
    return f
end

@everywhere function fg_gsfe_shared_1!(	storage_f::Union{Nothing, SharedArray{Float64,1}},
				                        storage_g::Union{Nothing, SharedArray{Float64,3}},
                                        M::Int64,
                                        invtE::Array{Array{Float64,2},1},
                                        Γ0::SharedArray{Float64,2},
                                        u1::SubArray{Float64,2},
                                        u2::SubArray{Float64,2},
                                        u3::SubArray{Float64,2},
                                        permutations::SharedArray{Int64,2},
                                        cache_misfit::SharedArray{Float64,1},
                                        cache_shifts::SharedArray{Float64,2},
                                        cache_gsfe_computation::SharedArray{Float64,2})
    idx = indexpids(storage_f)
    rge = myrange(idx, M)

    tmp = (storage_f != nothing ? 0.0 : nothing)
    if storage_g != nothing
        @inbounds for i=rge storage_g[:,i,1] .= 0.0 end
    end
@views @inbounds    tmp = gsfe_kernel(rge, tmp, storage_g[:,:,1], invtE[2],
                        Γ0, u1, u2, permutations[:,1],
                        cache_misfit, cache_shifts, cache_gsfe_computation )
@views @inbounds    tmp = gsfe_kernel(rge, tmp, storage_g[:,:,1], invtE[1],
                        Γ0, u1, u2, permutations[:,1],
                        cache_misfit, cache_shifts, cache_gsfe_computation )

    if storage_g != nothing
        @inbounds for i=rge storage_g[:,i,3] .= 0.0 end
    end
@views @inbounds    tmp = gsfe_kernel(rge, tmp, storage_g[:,:,3], invtE[2],
                        Γ0, u3, u2, permutations[:,2],
                        cache_misfit, cache_shifts, cache_gsfe_computation )
@views @inbounds    tmp = gsfe_kernel(rge, tmp, storage_g[:,:,3], invtE[3],
                        Γ0, u3, u2, permutations[:,2],
                        cache_misfit, cache_shifts, cache_gsfe_computation )

    @inbounds storage_f[idx] = 0.5*tmp
end


@everywhere function g_gsfe_shared_2!(storage::SharedArray{Float64,3},
                                            M::Int64,
                                            permutations::SharedArray{Int64,2})
    idx = indexpids(storage)
    rge = myrange(idx, M)

    @inbounds for i=rge
        storage[1,permutations[i,1],2] = -storage[1,i,1]
        storage[2,permutations[i,1],2] = -storage[2,i,1]
    end
end

@everywhere function g_gsfe_shared_3!(storage::SharedArray{Float64,3},
                                            M::Int64,
                                            permutations::SharedArray{Int64,2})
    idx = indexpids(storage)
    rge = myrange(idx, M)

    @inbounds for i=rge
        storage[1,permutations[i,2],2] -= storage[1,i,3]
        storage[2,permutations[i,2],2] -= storage[2,i,3]
    end
end

function EnergyGradient!(F::Union{Nothing, Float64},
                         G::Union{Nothing, Array{ComplexF64, 6}},
                         u::Array{ComplexF64, 6},
                         hull::Hull)
    N = hull.N
    hN = hull.hN
    # @timeit to "copy"
    copyto!(hull.cache_u, u)                       # Because FFTW destroys the array on input, even for out-of-place transform.

    @timeit to "g_elastic_shared" @sync begin
        for pw in workers()
            @async remotecall_wait(g_elastic_shared!, pw,
                        reshape(hull.cache_grad, (2, hN*N^3*3)),
                        hN*N^3*3, hull.Hess_Elastic, reshape(hull.cache_u, (2,hN*N^3*3)))
        end
    end
    @timeit to "f_elastic_shared" if F != nothing
        inc = 2*hN
        len = N^3*3
        F = real(
              BLAS.dotc(len*inc, pointer(u), 1, pointer(hull.cache_grad), 1)
        - 0.5*BLAS.dotc(len, pointer(u,1), inc, pointer(hull.cache_grad), inc)
        - 0.5*BLAS.dotc(len, pointer(u,2), inc, pointer(hull.cache_grad)+sizeof(ComplexF64), inc) )
    end

    @timeit to "backward fft" mul!(hull.cache_ur, hull.iplan, hull.cache_u)   # Multithreaded FFTW plan
    @timeit to "views" u1 = view(reshape(hull.cache_ur, (2,N^4,3)), :, :, 1)
    @timeit to "views" u2 = view(reshape(hull.cache_ur, (2,N^4,3)), :, :, 2)
    @timeit to "views" u3 = view(reshape(hull.cache_ur, (2,N^4,3)), :, :, 3) # displacement of each layer

    @timeit to "fg_gsfe_shared" @sync begin
        for pw in workers()
            @async remotecall_wait(fg_gsfe_shared_1!, pw,
                        F == nothing ? nothing : hull.cache_Misfits,
                        G == nothing ? nothing : hull.cache_grad_misfit,
                        N^4, hull.tl.invtE, hull.Γ0, u1, u2, u3,
                        hull.permutations, hull.cache_misfit, hull.cache_shifts, hull.cache_gsfe_computation)
        end
    end
    @timeit to "g_gsfe_shared" if G != nothing
        @sync begin
            for pw in workers()
                @async remotecall_wait(g_gsfe_shared_2!, pw,
                            hull.cache_grad_misfit, N^4, hull.permutations)
            end
        end
        @sync begin
            for pw in workers()
                @async remotecall_wait(g_gsfe_shared_3!, pw,
                            hull.cache_grad_misfit, N^4, hull.permutations)
            end
        end
    end
    if G != nothing
        @timeit to "forward fft" mul!(G,  hull.plan, reshape(hull.cache_grad_misfit, (2,N,N,N,N,3)))  # Multithreaded FFTW plan
        @timeit to "sum gradients" G .+= hull.cache_grad

#     # We need to enforce that the gradient has zero average,
#     # as the misfit energy is not exactly invariant under individual translations of the layers.
        G[:,1,1,1,1,:] .= 0.0
    end

    if F != nothing
        return (F + sum(hull.cache_Misfits)) / N^4
    end
    nothing
end

function Energy(u::Array{ComplexF64, 6}, hull::Hull)
    return EnergyGradient!(0.0, nothing, u, hull)
end

function Gradient!(grad::Array{ComplexF64, 6}, u::Array{ComplexF64, 6}, hull::Hull)
    EnergyGradient!(nothing, grad, u, hull)
end

function E_misfit(u::Array{Float64, 6}, hull::Hull)
    N = hull.N
    Γ0 = hull.Γ0
    invE = hull.tl.invtE
    misfit = 0
    u2 = reshape(u[:, :, :, :, :, 2], (2, N^4))

    for j = 1:2
        invE1 = invE[2]
        invE2 = invE[(j-1)*2+1]
        shifts1 = zeros(Float64, size(Γ0))
        shifts2 = zeros(Float64, size(Γ0))
        uj = reshape(u[:, :, :, :, :, (j-1)*2+1], (2, N^4))
        permutation = hull.permutations[:, j]
        shifts1[1,:] = Γ0[1,:] + invE1[1,1]*(u2[1, permutation] - uj[1,:]) + invE1[1,2]*(u2[2, permutation] - uj[2,:])
        shifts1[2,:] = Γ0[2,:] + invE1[2,1]*(u2[1, permutation] - uj[1,:]) + invE1[2,2]*(u2[2, permutation] - uj[2,:])

        shifts2[1,:] = Γ0[1,:] + invE2[1,1]*(u2[1, permutation] - uj[1,:]) + invE2[1,2]*(u2[2, permutation] - uj[2,:])
        shifts2[2,:] = Γ0[2,:] + invE2[2,1]*(u2[1, permutation] - uj[1,:]) + invE2[2,2]*(u2[2, permutation] - uj[2,:])
        misfit += 0.5 * sum(GSFE(shifts1)+GSFE(shifts2)) / N^4
    end

    return misfit
end

function E_elastic(u::Array{Float64, 6}, hull::Hull)
    N = hull.N
    hN = hull.hN
    Hess = hull.Hess_Elastic
    uk = hull.plan*u
    uk = reshape(uk,  (2,hN*N^3*3))

    graduH = zeros(ComplexF64, size(uk))

    for i = 1:hN*N^3*3
        graduH[1,i] = Hess[1,1,i] * uk[1,i] + Hess[2,1,i] * uk[2,i]
        graduH[2,i] = Hess[1,2,i] * uk[1,i] + Hess[2,2,i] * uk[2,i]
    end

    #graduH = hull.cache_grad
    # println(sum(graduH))
    inc = 2*hN
    len = N^3*3

    return real(
          BLAS.dotc(len*inc, pointer(uk), 1, pointer(graduH), 1)
    - 0.5*BLAS.dotc(len, pointer(uk,1), inc, pointer(graduH), inc)
    - 0.5*BLAS.dotc(len, pointer(uk,2), inc, pointer(graduH)+sizeof(ComplexF64), inc) ) / N^4
end


function Energy(u::Array{Float64, 6}, hull::Hull)
    return E_misfit(u, hull) + E_elastic(u, hull)
end

# Checked
function grad_misfit(u, hull)
    N = hull.N
    Γ0 = hull.Γ0
    invE = hull.tl.invtE
    misfit = 0
    u2 = reshape(u[:, :, :, :, :, 2], (2, N^4))
    grad = zeros(Float64, (2, N^4, 3))

    # GSFE gradient for layers 1 and 3
    for j = 1:2
        invE1 = invE[2]
        invE2 = invE[(j-1)*2+1]
        shifts1 = zeros(Float64, size(Γ0))
        shifts2 = zeros(Float64, size(Γ0))
        uj = reshape(u[:, :, :, :, :, (j-1)*2+1], (2, N^4))
        permutation = hull.permutations[:, j]
        shifts1[1,:] = Γ0[1,:] + invE1[1,1]*(u2[1, permutation] - uj[1,:]) + invE1[1,2]*(u2[2, permutation] - uj[2,:])
        shifts1[2,:] = Γ0[2,:] + invE1[2,1]*(u2[1, permutation] - uj[1,:]) + invE1[2,2]*(u2[2, permutation] - uj[2,:])

        gsfe_tmp = gradient_GSFE(shifts1)
        grad[1,:,2*j-1] -= 0.5 * (invE1[1,1] * gsfe_tmp[1,:]
                                + invE1[2,1] * gsfe_tmp[2,:])
        grad[2,:,2*j-1] -= 0.5 * (invE1[1,2] * gsfe_tmp[1,:]
                                + invE1[2,2] * gsfe_tmp[2,:])

        shifts2[1,:] = Γ0[1,:] + invE2[1,1]*(u2[1, permutation] - uj[1,:]) + invE2[1,2]*(u2[2, permutation] - uj[2,:])
        shifts2[2,:] = Γ0[2,:] + invE2[2,1]*(u2[1, permutation] - uj[1,:]) + invE2[2,2]*(u2[2, permutation] - uj[2,:])

        gsfe_tmp = gradient_GSFE(shifts2)
        grad[1,:,2*j-1] -= 0.5 * (invE2[1,1] * gsfe_tmp[1,:]
                                + invE2[2,1] * gsfe_tmp[2,:])
        grad[2,:,2*j-1] -= 0.5 * (invE2[1,2] * gsfe_tmp[1,:]
                                + invE2[2,2] * gsfe_tmp[2,:])
    end

    # permutation for the second layer
    for j = 1:2
        grad[j, hull.permutations[:,1],2] = -grad[j,:,1]
        grad[j, hull.permutations[:,2],2] -= grad[j,:,3]
    end

    grad = reshape(grad, (2,N,N,N,N,3))

    return grad_k = hull.plan*grad
end


function grad_elastic(u, hull)
    Hess = hull.Hess_Elastic
    N = hull.N
    hN = hull.hN
    uk = hull.plan*u
    uk = reshape(uk,  (2,hN*N^3*3))
    graduH = zeros(ComplexF64, size(uk))

    for i = 1:hN*N^3*3
        graduH[1,i] = Hess[1,1,i] * uk[1,i] + Hess[2,1,i] * uk[2,i]
        graduH[2,i] = Hess[1,2,i] * uk[1,i] + Hess[2,2,i] * uk[2,i]
    end

    return reshape(graduH, (2, hN, N, N, N, 3))
end

# input: u in real space (configuration), output energy gradient in k space
function EnergyGradient(u::Array{Float64, 6}, hull::Hull)
    return grad_misfit(u, hull) + grad_elastic(u, hull)
end


nothing
