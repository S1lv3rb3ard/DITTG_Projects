# One-time setup for example2.jl / example_vertical.jl.
#
#   julia install_julia_deps.jl
#
# Creates a local Julia project (Project.toml + Manifest.toml) in this folder
# so the relaxation code always gets the package versions it was written for.
# example2.jl uses the Optim 1.x API (Optim.only_fg!, LBFGS(P=...)), which was
# removed in Optim 2.x, so Optim is pinned to 1.x here.
#
# Afterwards run the relaxation with the project flag, e.g.
#   julia --project=. example2.jl -- 1.4 -2.8 48
using Pkg
Pkg.activate(@__DIR__)
Pkg.add([PackageSpec(name = "Optim", version = "1"),
         PackageSpec(name = "LineSearches"),
         PackageSpec(name = "ArgParse"),
         PackageSpec(name = "FFTW"),
         PackageSpec(name = "TimerOutputs"),
         PackageSpec(name = "Distributed"),
         PackageSpec(name = "SharedArrays"),
         PackageSpec(name = "LinearAlgebra"),
         PackageSpec(name = "Printf")])
Pkg.compat("Optim", "1")
Pkg.precompile()
Pkg.status()
println("Julia environment for RelaxedTTG is ready: use  julia --project=. example2.jl -- theta12 theta23 N")
