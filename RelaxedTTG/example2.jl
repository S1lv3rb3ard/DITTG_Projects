# Energy minimization
using Distributed
using TimerOutputs
using FFTW
using LinearAlgebra

using Optim
using LineSearches
using ArgParse
using Printf


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "theta12"
            help = "first twist angle theta12"
            required = true
            arg_type = Float64
        "theta23"
            help = "second twist angle theta23"
            required = true
            arg_type = Float64
        "N" 
            help = "discretization"
            required = true
            arg_type = Int64
end

    return parse_args(s)
end

parsed_args = parse_commandline()
println("Parsed args:")
for (arg,val) in parsed_args
    println("  $arg  =>  $val")
end

θ12_deg = parsed_args["theta12"]
θ23_deg = parsed_args["theta23"]

θ1 = deg2rad(θ12_deg);
θ2 = deg2rad(0.0);       # The second layer is fixed as a reference.
θ3 = deg2rad(θ23_deg);
N = parsed_args["N"]

global to = TimerOutput()

@timeit to "FFTW wisdom i/o" isfile("FFTWwisdom.jld") && FFTW.import_wisdom("FFTWwisdom.jld")
FFTW.set_num_threads(4)
BLAS.set_num_threads(4)

pcs = addprocs(4)
println("\n=============================================================\n        Using "
        * string(nworkers()) *  " workers... \n")

@timeit to "Trilayers.jl" include("Trilayers.jl")
@timeit to "GrapheneParameters.jl" @everywhere include("GrapheneParameters.jl")


@timeit to "Trilayer setup" tlg = Trilayer(l*E0, l*P0, θ1, θ3, K, G)
@timeit to "Hull setup" hull = Hull(tlg, N)
@timeit to "FFTW wisdom i/o" FFTW.export_wisdom("FFTWwisdom.jld")
hN = hull.hN;

f(u::Array{ComplexF64,6})                         = @timeit to "f(u)" Energy(u, hull)
g!( storage::Array{ComplexF64,6},
    u::Array{ComplexF64,6})                       = @timeit to "g!(storage, u)" Gradient!(storage, u, hull)
fg!(    F::Union{Nothing, Float64},
        G::Union{Nothing, Array{ComplexF64,6}},
        u::Array{ComplexF64,6})                   = @timeit to "fg!(F,G,u)" EnergyGradient!(F, G, u, hull)

@timeit to "Initialization" guess                   = zeros(ComplexF64, 2,hN,N,N,N,3)
@timeit to "Initialization" G = similar(guess)
# Precompilation
EnergyGradient!(0.0, G, guess, hull);
dot(G, guess);

# Run optimization algo
@timeit to "function setup" d = Optim.only_fg!(fg!)
@timeit to "method" method = LBFGS(; m = 5, P = hull.Precon_Elastic, scaleinvH0 = false)
@timeit to "options" options = Optim.Options(
            iterations = 10000, x_tol = 1e-4, f_tol = 1e-11, g_tol = 1e-8, allow_f_increases = true,
            show_trace = true, store_trace = false, show_every = 10 )
@timeit to "optimize" results = optimize(d, guess, method, options)
display(results)
display(to)

u = hull.iplan * results.minimizer
∇u = ∇(results.minimizer, hull)

outpath = string("$(@__DIR__)", "/data/");

write(string(outpath, @sprintf "triG_data_%.2f_%.2f_%d.jld" rad2deg(θ1) rad2deg(θ3) N), N, θ1, θ3, tlg.E, tlg.P, K, G)
write(string(outpath, @sprintf "triG_minimizer_%.2f_%.2f_%d.jld" rad2deg(θ1) rad2deg(θ3) N), u)
write(string(outpath, @sprintf "triG_gradient_%.2f_%.2f_%d.jld" rad2deg(θ1) rad2deg(θ3) N), ∇u[1,1], ∇u[1,2], ∇u[2,1], ∇u[2,2])

rmprocs(pcs)
nothing