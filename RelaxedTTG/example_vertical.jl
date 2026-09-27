# Corrugated, pressure-consistent trilayer relaxation.
#
#   julia example_vertical.jl <theta12> <theta23> <N> [pressure_GPa] [kappa_eV]
#
# Writes to data/:
#   triG_vz_data_<t12>_<t23>_<N>_P<P>.jld       same header as triG_data (N, theta1,
#                                               theta3, E, P, K, G)
#   triG_vz_minimizer_<t12>_<t23>_<N>_P<P>.jld  u, same layout as triG_minimizer
#   triG_vz_vertical_<t12>_<t23>_<N>_P<P>.jld   pressure, kappa, dbar12, dbar23, w
# (raw little-endian Base.write output, like example2.jl). A negative angle may
# need `--` before the positional arguments.

using Optim
using Printf

include("TrilayersVertical.jl")

args = ARGS
length(args) >= 3 || error("usage: julia example_vertical.jl theta12 theta23 N [P_GPa] [kappa]")
θ12 = parse(Float64, args[1])
θ23 = parse(Float64, args[2])
N = parse(Int, args[3])
pressure = length(args) >= 4 ? parse(Float64, args[4]) : 0.0
κ = length(args) >= 5 ? parse(Float64, args[5]) : 1.4

h = VerticalHull(θ12, θ23, N; κ = κ, pressure = pressure)
ε0 = -log(1 + pressure / PRESSURE_A) / PRESSURE_B
M = N^4
x0 = zeros(2 * M * 3 + M * 3 + 2)
x0[end-1:end] .= D0 * (1 + ε0)

fg!(F, G, x) = energy_gradient!(G, x, h)
result = optimize(Optim.only_fg!(fg!), x0,
                  LBFGS(; m = 30),
                  Optim.Options(iterations = 30000, g_tol = 1e-11, f_tol = 1e-16,
                                show_trace = true, show_every = 50))
display(result)
x = Optim.minimizer(result)
u = reshape(x[1:2*M*3], (2, N, N, N, N, 3))
w = reshape(x[2*M*3+1:2*M*3+M*3], (N, N, N, N, 3))
dbar = x[end-1:end]
@printf("dbar12 = %.4f A, dbar23 = %.4f A, max|w| = %.4f A, max|u| = %.4f A\n",
        dbar[1], dbar[2], maximum(abs, w), maximum(abs, u))

outpath = string(@__DIR__, "/data/")
mkpath(outpath)
tag = @sprintf("%.2f_%.2f_%d_P%.2f", θ12, θ23, N, pressure)
write(string(outpath, "triG_vz_data_", tag, ".jld"), N, h.θ[1], h.θ[3], h.E, h.P, K, G)
write(string(outpath, "triG_vz_minimizer_", tag, ".jld"), u)
write(string(outpath, "triG_vz_vertical_", tag, ".jld"), pressure, κ, dbar[1], dbar[2], w)
