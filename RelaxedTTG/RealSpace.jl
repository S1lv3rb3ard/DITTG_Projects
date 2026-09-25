## Given displacement discretized on a trilayer Hull,
# provide functions to compute relaxed and unrelaxed
# real space positions as well as resulting diffraction patterns.

include("Trilayers.jl")
# include("WSe2_0_Parameters.jl")
using Interpolations
using Printf

function Read(N, θ1, θ3, outpath, system)
  streamData = open(string( outpath, system,  "_data_", @sprintf "%.2f_%.2f_%d.jld" θ1 θ3 N))
  streamMinimizer = open(string( outpath, system, "_minimizer_", @sprintf "%.2f_%.2f_%d.jld" θ1 θ3 N))
  streamGradient = open(string( outpath, system, "_gradient_", @sprintf "%.2f_%.2f_%d.jld" θ1 θ3 N))


  N = read(streamData, Int64)
  θ1 = read(streamData, Float64)
  θ3 = read(streamData, Float64)
  E = zeros(Float64,2,2); read!(streamData, E)
  P = zeros(Float64,2); read!(streamData, P)
  #K = read(streamData, Float64)
  #G = read(streamData, Float64)

  u = zeros(Float64,2,N,N,N,N,3); read!(streamMinimizer, u)
  ux1 = zeros(Float64,N,N,N,N,3); read!(streamGradient, ux1)
  uy1 = zeros(Float64,N,N,N,N,3); read!(streamGradient, uy1)
  # ∇u = [[ux1] [uy1]]
  ux2 = zeros(Float64,N,N,N,N,3); read!(streamGradient, ux2)
  uy2 = zeros(Float64,N,N,N,N,3); read!(streamGradient, uy2)
  ∇u = [[ux1] [uy1];
        [ux2] [uy2]]
  close(streamData)
  close(streamMinimizer)
  close(streamGradient)
  println(string("G = ", G))
  println(string("K = ", K))
  tl = Trilayer(E, P, θ1, θ3, K, G);
  return tl, u, ∇u
end

mutable struct Configuration
    γ::Array{Float64,2}

    function Configuration(γ1::Array{Float64,1}, γ2::Array{Float64,1}, γ3::Array{Float64,1}, tl::Trilayer)
        γ = hcat(   tl.tE[1]*mod.(inv(tl.tE[1])*γ1, 1.),
                    tl.tE[2]*mod.(inv(tl.tE[2])*γ2, 1.),
                    tl.tE[3]*mod.(inv(tl.tE[3])*γ3, 1.))
        new(γ)
    end
end

mutable struct Displacement
    itp::Array{Interpolations.Extrapolation,2}     # Array of interpolants for the modulation functions

    function Displacement(u::Array{Float64,6}, N::Int)
        itp = Array{Interpolations.Extrapolation, 2}(undef,2,3)

        for i in 1:3
            itp[1,i] = extrapolate(interpolate(u[1,:,:,:,:,i], BSpline(Quadratic(Line())), OnGrid()), Periodic())
            itp[2,i] = extrapolate(interpolate(u[2,:,:,:,:,i], BSpline(Quadratic(Line())), OnGrid()), Periodic())
        end
        new(itp)
    end
end


# Generate the unrelaxed positions of atoms in a given radius for a given configuration

function create_positions(R::Float64, ω::Configuration, tl::Trilayer)
    n = floor(Int64, R/norm(inv(tl.E)))+1
    I = (-n:n) .* ones(Int64, (1,2*n+1))
    J = ones(Int64, (2*n+1,1)) .* (-n:n)'

    X1 = ω.γ[:,1] .+ tl.tE[1] * [I[:]'; J[:]'];
    X2 = ω.γ[:,2] .+ tl.tE[2] * [I[:]'; J[:]'];
    X3 = ω.γ[:,3] .+ tl.tE[3] * [I[:]'; J[:]'];

    X1 = X1[:, hypot.(X1[1,:],X1[2,:]) .< R]
    X2 = X2[:, hypot.(X2[1,:],X2[2,:]) .< R]
    X3 = X3[:, hypot.(X3[1,:],X3[2,:]) .< R]
    return [[X1]; [X2]; [X3]]
end

# Relax the position of atoms in both layers according to the provided displacement.
# Positions is an array of arrays of unrelaxed positions indexed by the layer index (1,2),
#    where the first dimension corresponds to directions x and y,
#          the second dimension to point index.
function displace!(Positions::Array{Array{Float64, 2}, 1}, ϕ::Displacement, ω::Configuration, N::Int, tl::Trilayer)
    # Layer 1
    CartPos12 = N .* (inv(tl.tE[2]) * (ω.γ[:,2] .- Positions[1])).+1
    CartPos13 = N .* (inv(tl.tE[3]) * (ω.γ[:,3] .- Positions[1])).+1
    for i=1:2, ind=1:size(Positions[1], 2)
        Positions[1][i,ind] += ϕ.itp[i,1](CartPos12[1,ind], CartPos12[2,ind], CartPos13[1,ind], CartPos13[2,ind])
    end
    # Layer 2
    CartPos21 = N .* (inv(tl.tE[1]) * (ω.γ[:,1] .- Positions[2])).+1
    CartPos23 = N .* (inv(tl.tE[3]) * (ω.γ[:,3] .- Positions[2])).+1
    for i=1:2, ind=1:size(Positions[2], 2)
        Positions[2][i,ind] += ϕ.itp[i,2](CartPos21[1,ind], CartPos21[2,ind], CartPos23[1,ind], CartPos23[2,ind])
    end
    # Layer 3
    CartPos31 = N .* (inv(tl.tE[1]) * (ω.γ[:,1] .- Positions[3])).+1
    CartPos32 = N .* (inv(tl.tE[2]) * (ω.γ[:,2] .- Positions[3])).+1
    for i=1:2, ind=1:size(Positions[3], 2)
        Positions[3][i,ind] += ϕ.itp[i,3](CartPos32[1,ind], CartPos32[2,ind], CartPos31[1,ind], CartPos31[2,ind])
    end
    nothing
end


function fourier_zoom(Positions::Array{Array{Float64,2},1},
                        sigma::Float64,
                        K::Array{Float64,1},
                        L::Float64,
                        N::Int64 = 64.,
                        m::Float64 = 1.)
    grid = (-N:N)

    filtereddata = zeros(Complex128, 2*N+1, 2*N+1)
    X = [0.;0.]
    for j=1:length(Positions)
        for index = 1:size(Positions[j],2)
            X[1] = Positions[j][1,index]/L;
            X[2] = Positions[j][2,index]/L;
            ϕ = exp(-2im*pi*L*dot(K,X) -.5*(L/sigma)^2*dot(X,X))
            for k = max(1, N+1 + floor(Int64, X[1]) - 6):min(2*N, N+1 + ceil(Int64, X[1]) + 6),
                l = max(1, N+1 + floor(Int64, X[2]) - 6):min(2*N, N+1 + ceil(Int64, X[2]) + 6)
                filtereddata[k,l] += ϕ * exp(-.5*((grid[k] - X[1])^2  + (grid[l] - X[2])^2))
            end
        end
    end

    filter = zeros(2*N+1)
    for k=max(1, N-6):min(2*N+1, N+6)
        filter[k] = exp(-.5*grid[k]^2)
    end
    intensity = abs2.(fftshift(fft(filter)))

    S = sum([sum(exp.( (-.5/sigma^2) .* sum(Pos.^2, 1))) for Pos in Positions])
    farfieldintensity = abs2.(fftshift(fft(filtereddata))) ./ intensity ./ intensity' ./ S


    zoom = 1+floor(Int64, 2*N/3):ceil(Int64, 4*N/3)
    Kx = K[1] .+ (m/(L*(2*N+1))) .* grid[zoom]' .* ones(zoom)
    Ky = K[2] .+ (m/(L*(2*N+1))) .* ones(zoom') .* grid[zoom]
    I = farfieldintensity[zoom,zoom]'

    return Kx, Ky, I
end


mutable struct Fields
    itp::Array{Interpolations.Extrapolation}     # Array of interpolants for the modulation functions
    shape
    num

    function Fields(g::Array{Array{Float64,5}}, N::Int)
        num = length(g)
        shape = size(g)
        itp = Array{Interpolations.Extrapolation}(undef, shape...,3)

        for j in 1:num, i in 1:3
            itp[j+num*(i-1)] = extrapolate(interpolate(g[j][:,:,:,:,i], BSpline(Quadratic(Periodic(OnCell())))), Periodic())
        end
        new(itp,shape,num)
    end
end

function interpolateFields(Positions::Array{Float64, 2}, ϕ::Fields, ω::Configuration, H::Int, tl::Trilayer)
    g = Array{Array{Float64,1}}(undef, ϕ.shape..., 3)
    for i in 1:3, n in 1:ϕ.num
        g[n+ϕ.num*(i-1)] = zeros(Float64, size(Positions,2))
    end

    # Layer 1
    CartPos12 = N .* (inv(tl.tE[2]) * ω.γ[:,2] .- inv(tl.tE[1]) * ω.γ[:,1]
                .+ (inv(tl.tE[1]) - inv(tl.tE[2])) * Positions).+1
    CartPos13 = N .* (inv(tl.tE[3]) * ω.γ[:,3] .- inv(tl.tE[1]) * ω.γ[:,1]
                .+ (inv(tl.tE[1]) - inv(tl.tE[3])) * Positions).+1
    for j in 1:ϕ.num, ind in 1:size(Positions, 2)
        g[j][ind] = ϕ.itp[j](CartPos12[1,ind], CartPos12[2,ind], CartPos13[1,ind], CartPos13[2,ind])
    end
    # Layer 2
    CartPos21 = N .* (inv(tl.tE[1]) * ω.γ[:,1] .- inv(tl.tE[2]) * ω.γ[:,2]
                .+ (inv(tl.tE[2]) - inv(tl.tE[1])) * Positions).+1
    CartPos23 = N .* (inv(tl.tE[3]) * ω.γ[:,3] .- inv(tl.tE[2]) * ω.γ[:,2]
                .+ (inv(tl.tE[2]) - inv(tl.tE[3])) * Positions).+1
    for j in 1:ϕ.num, ind in 1:size(Positions, 2)
        g[j+ϕ.num][ind] += ϕ.itp[j+ϕ.num](CartPos21[1,ind], CartPos21[2,ind], CartPos23[1,ind], CartPos23[2,ind])
    end
    # Layer 3
    CartPos31 = N .* (inv(tl.tE[1]) * ω.γ[:,1] .- inv(tl.tE[3]) * ω.γ[:,3]
                .+ (inv(tl.tE[3]) - inv(tl.tE[1])) * Positions).+1
    CartPos32 = N .* (inv(tl.tE[2]) * ω.γ[:,2] .- inv(tl.tE[3]) * ω.γ[:,3]
                .+ (inv(tl.tE[3]) - inv(tl.tE[2])) * Positions).+1
    for j in 1:ϕ.num, ind in 1:size(Positions, 2)
        g[j+2*ϕ.num][ind] += ϕ.itp[j+2*ϕ.num](CartPos32[1,ind], CartPos32[2,ind], CartPos31[1,ind], CartPos31[2,ind])
    end

    return g
end

#Reduce framework without allocation - Fabian Gans
#Define abstract Array which consumes values on setindex!
mutable struct ArrayReducer{T,N,F} <: AbstractArray{T,N}
    v::T
    size::NTuple{N,Int}
    op::F

    function ArrayReducer(size::NTuple{N},v0::T,op) where {T,N}
        return new{T,N,typeof(op)}(v0,size,op)
    end
end
ArrayReducer{T,N}(size::NTuple{N},v0::T,op) where {T,N} =ArrayReducer{T,N,typeof(op)}(v0,size,op)
Base.setindex!(x::ArrayReducer,v,i...) = x.v=x.op(x.v,v)
Base.size(a::ArrayReducer)=a.size
get_stop(a::Base.OneTo)=a.stop

#Define the function
function broadcast_reduce(f,op,v0,A,Bs...)
    shape = Base.Broadcast.combine_axes(A,Bs...)
    reducer = ArrayReducer(get_stop.(shape),v0,op)
    broadcast!(f,reducer,A,Bs...)
    reducer.v
end

function fft_coeffs(X::Array{Float64,1}, N::Int)
    hN = div(N,2)+1
    c1 = ones(ComplexF64, hN)
    c1[2:hN] .= 2.0.*exp.((2im*pi).*(1:hN-1)*X[1])
    c2 = exp.((2im*pi).*[0:hN-1;hN-N:-1]*X[2])
    c3 = exp.((2im*pi).*[0:hN-1;hN-N:-1]*X[3])
    c4 = exp.((2im*pi).*[0:hN-1;hN-N:-1]*X[4])

    if mod(N,2) == 0 # Adjust Nyquist frequency coefficient
        c1[hN] *= 0.5
    end

    return reshape(c1./N^2, (hN,1,1,1)), reshape(c2, (1,N,1,1)),
           reshape(c3, (1,1,N,1)) , reshape(c4, (1,1,1,N))
end


function evaluateFields(Positions::Array{Float64, 2}, ϕ::Array{Array{ComplexF64,5}}, ω::Configuration, N::Int, tl::Trilayer)
    hN = div(N,2)+1
    num = length(ϕ)
    g = Array{Array{Float64,1}}(undef, size(ϕ)..., 3)
    for i in 1:3, n in 1:num
        g[n+num*(i-1)] = zeros(Float64, size(Positions,2))
    end

    # Layer 1
    CartPos = cat(
                inv(tl.tE[2]) * ω.γ[:,2] .- inv(tl.tE[1]) * ω.γ[:,1]
                .+ (inv(tl.tE[1]) - inv(tl.tE[2])) * Positions,
                inv(tl.tE[3]) * ω.γ[:,3] .- inv(tl.tE[1]) * ω.γ[:,1]
                .+ (inv(tl.tE[1]) - inv(tl.tE[3])) * Positions,
                    dims = 1)
    for ind in 1:size(Positions, 2)
        c1, c2, c3, c4 = fft_coeffs(CartPos[:,ind], N)
        for n in 1:num
            g[n][ind] = real(broadcast_reduce(*,+,0.0im,c1,c2,c3,c4,view(ϕ[n],:,:,:,:,1)))
        end
    end

    # Layer 2
    CartPos = cat(
                inv(tl.tE[1]) * ω.γ[:,1] .- inv(tl.tE[2]) * ω.γ[:,2]
                .+ (inv(tl.tE[2]) - inv(tl.tE[1])) * Positions,
                inv(tl.tE[3]) * ω.γ[:,3] .- inv(tl.tE[2]) * ω.γ[:,2]
                .+ (inv(tl.tE[2]) - inv(tl.tE[3])) * Positions,
                    dims = 1)
    for ind in 1:size(Positions, 2)
        c1, c2, c3, c4 = fft_coeffs(CartPos[:,ind], N)
        for n in 1:num
            g[n+num][ind] = real(broadcast_reduce(*,+,0.0im,c1,c2,c3,c4,view(ϕ[n],:,:,:,:,2)))
        end
    end

    # Layer 3
    CartPos = cat(
                inv(tl.tE[2]) * ω.γ[:,2] .- inv(tl.tE[3]) * ω.γ[:,3]
                .+ (inv(tl.tE[3]) - inv(tl.tE[2])) * Positions,
                inv(tl.tE[1]) * ω.γ[:,1] .- inv(tl.tE[3]) * ω.γ[:,3]
                .+ (inv(tl.tE[3]) - inv(tl.tE[1])) * Positions,
                    dims = 1)
    for ind in 1:size(Positions, 2)
        c1, c2, c3, c4 = fft_coeffs(CartPos[:,ind], N)
        for n in 1:num
            g[n+2*num][ind] = real(broadcast_reduce(*,+,0.0im,c1,c2,c3,c4,view(ϕ[n],:,:,:,:,3)))
        end
    end

    return g
end

function evaluate(P::Array{Int,1}, ϕ::Array{ComplexF64,4}, N::Int)
    X = (1.0/N) .* (P.-1)
    c1, c2, c3, c4 = fft_coeffs(X, N)
    return real(broadcast_reduce(*,+,0.0im,c1,c2,c3,c4,ϕ))
end

nothing
