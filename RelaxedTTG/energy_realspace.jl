using LinearAlgebra, ArgParse

include("RealSpace.jl")
include("Trilayers.jl")

system = "triG"
include("GrapheneParameters.jl")


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
        "n"
            help = "output discretization"
            required = true 
            arg_type = Int64
        "scale"
            help = "scaling factor in the unit of moiré of moiré scale" 
            required = true 
            arg_type = Float64
        "path"
            help = "data file path"
            required = true
            arg_type = String
end

    return parse_args(s)

end

parsed_args = parse_commandline()
println("Parsed args:")
for (arg,val) in parsed_args
    println("  $arg  =>  $val")
end

θ1 = parsed_args["theta12"]
θ3 = parsed_args["theta23"]
N = parsed_args["N"]
n = parsed_args["n"]
scale = parsed_args["scale"]

outpath = string("$(@__DIR__)", parsed_args["path"])

println(θ1)
println(θ3)
println(N)

tlg, u, ∇u = Read(N, θ1, θ3, outpath, system)
hull = Hull(tlg, N)
Γ0 = Array(hull.Γ0)
permutations = Array(hull.permutations)
invE = hull.tl.invtE

u_tmp = reshape(u, (2, N^4, 3))

# Calculate the total relaxed energy 
Ee = E_elastic(u, hull)
Em = E_misfit(u, hull)
println( "Total elastic energy (eV/unit cell): ", Ee )
println( "Total misfti energy (eV/unit cell): ", Em )

# define function to calculate the moiré length 
function moireh_calc(A0, θ, δ, m, n)
    rot = [cos(θ) -sin(θ); sin(θ) cos(θ)] # ccw rotation
    A1 = A0
    A2 = (1+δ)*inv(rot)*A0

    g1 = 2*pi*inv(A1)'
    g2 = 2*pi*inv(A2)'

    G = m*g1 - n*g2
    Am = inv(G)' * (2*pi)
    return Am
end


function calc_moire_dom_tri(q12, q23, A0)
    grid_search = 8

    A12 = moireh_calc(A0, q12, 0, 1, 1)
    A23 = moireh_calc(A0, q23, 0, 1, 1)

    if q12*q23 > 0 
        rot180 = [cos(π) -sin(π); sin(π) cos(π)];
        A12 = rot180*A12;
    end 
    
    G12 = 2π*inv(A12)'
    G23 = 2π*inv(A23)'

    a12 = -A12[:,1]
    a23 = A23[:,1]

    g12 = -G12[:,1]
    g23 = G23[:,1]

    # bilayer moire length
    ml12 = norm(a12)
    ml23 = norm(a23)

    # angle between the bilayer moires 
    theta = acos(dot(a12, a23)/(ml12*ml23))

    # find the dominant harmonic 
    g_vec_norm = zeros(grid_search, grid_search)
    for i1 = 1:grid_search 
        for i2 = 1:grid_search 
            g_vec = i1*g12 - i2*g23 
            g_vec_norm[i1,i2] = norm(g_vec)
        end 
    end 
    
    idx = findall(g_vec_norm .== minimum(g_vec_norm[:]))[1]

    m_dom = idx[1]
    n_dom = idx[2]
    
    delta = ml12/ml23*m_dom/n_dom-1
    Am_dom = moireh_calc(-A12, theta, ml23/ml12-1, m_dom, n_dom)
    return Am_dom
end 


# relaxed energy
# L1 + L2
shifts = zeros(Float64, (2,N^4))
uj = u_tmp[:,:,1]
u2 = u_tmp[:,:,2]
permutation = permutations[:,1]
for i = 1:2
    if i == 1
        global misfit12 = zeros(Float64, (N^4))
    end

    invEi = invE[i]
    shifts[1,:] = Γ0[1,:] + invEi[1,1]*(u2[1, permutation] - uj[1,:]) +
                                  invEi[1,2]*(u2[2, permutation] - uj[2,:])
    shifts[2,:] = Γ0[2,:] + invEi[2,1]*(u2[1, permutation] - uj[1,:]) +
                                  invEi[2,2]*(u2[2, permutation] - uj[2,:])
    misfit12 = misfit12 + 0.5*GSFE(shifts)
end

# L2 + L3
uj = u_tmp[:,:,2]
u2 = u_tmp[:,:,3]
permutation = permutations[:,2]
for i = 1:2
    if i == 1
        global misfit23 = zeros(Float64, (N^4))
    end
    invEi = invE[i+1]
    shifts[1,:] = Γ0[1,:] + invEi[1,1]*(u2[1, permutation] - uj[1,:]) +
                                  invEi[1,2]*(u2[2, permutation] - uj[2,:])
    shifts[2,:] = Γ0[2,:] + invEi[2,1]*(u2[1, permutation] - uj[1,:]) +
                                  invEi[2,2]*(u2[2, permutation] - uj[2,:])
    misfit23 += 0.5*GSFE(shifts)
end

# unrelaxed energy
# L1 + L2
shifts = zeros(Float64, (2,N^4))
uj = u_tmp[:,:,1]
u2 = u_tmp[:,:,2]
permutation = permutations[:,1]
for i = 1:2
    if i == 1
        global misfit12_unrelaxed = zeros(Float64, (N^4))
    end

    invEi = invE[i]
    shifts[1,:] = Γ0[1,:] 
    shifts[2,:] = Γ0[2,:] 
    misfit12_unrelaxed = misfit12_unrelaxed + 0.5*GSFE(shifts)
end

# L2 + L3
uj = u_tmp[:,:,2]
u2 = u_tmp[:,:,3]
permutation = permutations[:,2]
for i = 1:2
    if i == 1
        global misfit23_unrelaxed = zeros(Float64, (N^4))
    end
    invEi = invE[i+1]
    shifts[1,:] = Γ0[1,:] 
    shifts[2,:] = Γ0[2,:] 
    misfit23_unrelaxed += 0.5*GSFE(shifts)
end

Γx = reshape(Γ0[1, :], (N, N, N, N))
Γy = reshape(Γ0[2, :], (N, N, N, N))
b = zeros(Float64, (2, N, N))
b[1, :, :] = Γx[:, :, 1, 1]
b[2, :, :] = Γy[:, :, 1, 1]
b = reshape(b, (2, N^2))
b = tlg.E * b
bx = reshape(b[1, :], (N, N))
by = reshape(b[2, :], (N, N))

misfit12 = reshape(misfit12, (N, N, N, N))
misfit23 = reshape(misfit23, (N, N, N, N))

misfit12 = reshape(misfit12, (N,N,N,N,1))
misfit23 = reshape(misfit23, (N,N,N,N,1))
misfit = zeros((N,N,N,N,3))
misfit[:,:,:,:,1]=misfit12
misfit[:,:,:,:,3]=misfit23
misfit12 = reshape(misfit12, (N, N, N, N))
misfit23 = reshape(misfit23, (N, N, N, N))
itp_misfit = Fields([misfit], N)

# unrelaxed 
misfit12_unrelaxed = reshape(misfit12_unrelaxed, (N,N,N,N,1))
misfit23_unrelaxed = reshape(misfit23_unrelaxed, (N,N,N,N,1))
misfit_unrelaxed = zeros((N,N,N,N,3))
misfit_unrelaxed[:,:,:,:,1]=misfit12_unrelaxed
misfit_unrelaxed[:,:,:,:,3]=misfit23_unrelaxed
misfit12_unrelaxed = reshape(misfit12_unrelaxed, (N,N,N,N))
misfit23_unrelaxed = reshape(misfit23_unrelaxed, (N,N,N,N))
itp_misfit_unrelaxed = Fields([misfit_unrelaxed], N)

itp_u = Fields([u[1,:,:,:,:,:], u[2,:,:,:,:,:]], N)
itp_∇u = Fields(∇u, N)

# E = ( abs(tlg.θ[1]) < abs(tlg.θ[3]) ?
#             inv(inv(tlg.tE[1]) - inv(tlg.tE[2])) :
#             inv(inv(tlg.tE[3]) - inv(tlg.tE[2])) )

# Lx = maximum(abs.(E[1,:]))
# Ly = maximum(abs.(E[2,:]))
ω = Configuration([0.0;0.0], [0.0;0.0], [0.0;0.0], tlg)

A = [1 sqrt(3)/2; 0 1/2] * sqrt(3) * l 
A_sc = calc_moire_dom_tri(deg2rad(θ1), deg2rad(θ3), A )
a_sc = norm(A_sc[:,1])

a_sc = a_sc * scale

X = range(-a_sc, a_sc, length=4*n)
Y = range(-a_sc, a_sc*1.5, length=5*n)
Grid = cat( reshape(X, (1,4*n,1)) .* ones(1,1,5*n),
            ones(1,4*n,1) .* reshape(Y, (1,1,5*n)), dims=1)
Grid = reshape(Grid, (2, 20*n^2))
# @time exactU = evaluateFields(Grid, ϕ, ω, N, tlg)
@time U = interpolateFields(Grid, itp_u, ω, N, tlg)
@time ∇U = interpolateFields(Grid, itp_∇u, ω, N, tlg)
@time E_interp = interpolateFields(Grid, itp_misfit, ω, N, tlg)
@time E_interp_unrelaxed = interpolateFields(Grid, itp_misfit_unrelaxed, ω, N, tlg)


# close("all")
xarr = X .* ones(1,5*n)
yarr = ones(4*n,1) .* reshape(Y, (1,5*n))

E_interp12 = reshape(E_interp[1], size(xarr))
E_interp23 = reshape(E_interp[3], size(xarr))
E12 = reshape(E_interp_unrelaxed[1], size(xarr))
E23 = reshape(E_interp_unrelaxed[3], size(xarr))

xarr = xarr[:]
yarr = yarr[:]
E_interp12 = E_interp12[:]
E_interp23 = E_interp23[:]
E12 = E12[:]
E23 = E23[:]

grad1 = ∇U[2,1,1] .- ∇U[1,2,1]
grad2 = ∇U[2,1,2] .- ∇U[1,2,2]
grad3 = ∇U[2,1,3] .- ∇U[1,2,3]

open(string(outpath, system, "_q12_", round(θ1; digits=2), "deg_q23_",
    round(θ3; digits=2), "deg_N_", N, "_scale_", scale, "_n_", n, "_disp_energy.txt"), "w") do f
    write(f, "x,y,u1x,u1y,u2x,u2y,u3x,u3y,misfit12,misfit23,curl1,curl2,curl3,misfit12_urlx,misfit23_urlx \n")
    for i in 1:length(E_interp12)
        rx = xarr[i]
        ry = yarr[i]
        r1 = U[1, 1][i]
        r2 = U[2, 1][i]
        r3 = U[1, 2][i]
        r4 = U[2, 2][i]
        r5 = U[1, 3][i]
        r6 = U[2, 3][i]
        r7 = E_interp12[i]
        r8 = E_interp23[i]
        r9 = grad1[i] 
        r10 = grad2[i]
        r11 = grad3[i]
        r12 = E12[i]
        r13 = E23[i]
        write(f, "$rx, $ry, $r1, $r2, $r3, $r4, $r5, $r6, $r7, $r8, $r9, $r10, $r11, $r12, $r13\n")
    end
end





