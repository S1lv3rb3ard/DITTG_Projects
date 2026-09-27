# Pressure-consistent relaxation of a twisted trilayer with out-of-plane corrugation.
#
# Extends the configuration-space model of Trilayers.jl (Cazeaux, Zhu et al.) with
#   * out-of-plane fields w_j on the same hulls (zero mean) and the mean interlayer
#     spacings dbar12, dbar23,
#   * the bending energy (1/2) kappa Omega <(Lap_r w_j)^2>,
#   * a distance-dependent stacking energy
#         Phi(s, d) = c0(d) + R(d) * GSFE(s),
#     where GSFE is the zero-pressure functional of GrapheneParameters.jl,
#     R(d) = exp(poly(eps)) is the distance scaling of the AA-AB energy of the refined
#     Kolmogorov-Crespi potential (Ouyang et al., Nano Lett. 18, 6009 (2018)), and c0
#     is the registry-independent binding energy,
#   * the enthalpy P Omega (dbar12 + dbar23) at uniaxial pressure P.
# eps = d/d0 - 1 with d0 = 3.35 A (negative under compression). The default binding
# energy integrates the pressure law P = A (exp(-B eps) - 1) of Carr et al.,
# PRB 98, 085144 (2018), so the pressure scale matches the interlayer hopping fits.
#
# Energies per unit cell (eV), lengths in Angstrom, pressure in GPa. The model and its
# gradient were verified against finite differences in a Python port of this file.
#
# Layout: u[c, s, t, v, w, layer] exactly as Trilayers.jl; w[s, t, v, w, layer].
# Layer j stores its configuration pairs relative to layers (2,3), (1,3), (2,1).

using FFTW
using LinearAlgebra

include("GrapheneParameters.jl")    # l, E0, P0, K, G, GSFE, gradient_GSFE

const D0 = 3.35
const GPaA3 = 6.241509e-3            # eV per GPa*A^3
const PRESSURE_A = 5.73              # GPa
const PRESSURE_B = 9.54
const LOGR = [7.25138301e-02, 3.83435782e-02, -1.07401220e+01, 0.0]   # eps^3 .. eps^0

evalpoly_desc(p, x) = foldl((acc, c) -> acc .* x .+ c, p; init = zero(x))
derivpoly_desc(p) = [p[i] * (length(p) - i) for i in 1:length(p)-1]

struct VerticalHull
    N::Int
    θ::Vector{Float64}
    E::Matrix{Float64}
    P::Vector{Float64}
    invtE::Vector{Matrix{Float64}}     # 2π inv(tE_j)
    Γ0::Array{Float64,5}               # (2,N,N,N,N) hull disregistry in radians
    perm1::Array{CartesianIndex{4},4}  # layer-1 point -> layer-2 point
    perm2::Array{CartesianIndex{4},4}  # layer-3 point -> layer-2 point
    tK::Array{Float64,6}               # (2,N,N,N,N,3) real-space gradient multipliers
    C4::Array{Float64,4}
    Ω::Float64
    κ::Float64
    pressure::Float64
end

function jrot(t)
    [cos(t) sin(t); -sin(t) cos(t)]
end

function VerticalHull(θ12deg, θ23deg, N; κ = 1.4, pressure = 0.0)
    E = l * E0
    P = l * P0
    θ = [deg2rad(θ12deg), 0.0, deg2rad(θ23deg)]
    tE = [jrot(t) * E for t in θ]
    inv_tE = inv.(tE)
    invtE = [2π * M for M in inv_tE]
    Γ0 = zeros(2, N, N, N, N)
    perm1 = Array{CartesianIndex{4},4}(undef, N, N, N, N)
    perm2 = Array{CartesianIndex{4},4}(undef, N, N, N, N)
    m(a) = mod(a, N) + 1
    for w in 0:N-1, v in 0:N-1, t in 0:N-1, s in 0:N-1
        Γ0[1, s+1, t+1, v+1, w+1] = 2π * s / N
        Γ0[2, s+1, t+1, v+1, w+1] = 2π * t / N
        perm1[s+1, t+1, v+1, w+1] = CartesianIndex(m(-s), m(-t), m(v - s), m(w - t))
        perm2[s+1, t+1, v+1, w+1] = CartesianIndex(m(v - s), m(w - t), m(-s), m(-t))
    end
    T1 = [(inv_tE[1] - inv_tE[2])', (inv_tE[2] - inv_tE[1])', (inv_tE[3] - inv_tE[2])']
    T2 = [(inv_tE[1] - inv_tE[3])', (inv_tE[2] - inv_tE[3])', (inv_tE[3] - inv_tE[1])']
    k = 2π .* [0:(N-1)÷2; -(N÷2):-1]         # integer frequencies (FFT order) times 2π
    tK = zeros(2, N, N, N, N, 3)
    for L in 1:3, w in 1:N, v in 1:N, t in 1:N, s in 1:N, i in 1:2
        tK[i, s, t, v, w, L] = T1[L][i, 1] * k[s] + T1[L][i, 2] * k[t] +
                               T2[L][i, 1] * k[v] + T2[L][i, 2] * k[w]
    end
    Cmat = [K+G 0 0 K-G; 0 G G 0; 0 G G 0; K-G 0 0 K+G]
    VerticalHull(N, θ, E, P, invtE, Γ0, perm1, perm2, tK,
                 reshape(Cmat, (2, 2, 2, 2)), abs(det(E)), κ, pressure)
end

# Registry-independent binding energy c0(d) and its derivative (eV per cell, eV/A).
function binding(h::VerticalHull, d)
    e = d ./ D0 .- 1
    K0 = h.Ω * D0 * PRESSURE_A * GPaA3
    c0 = K0 .* ((exp.(-PRESSURE_B .* e) .- 1) ./ PRESSURE_B .+ e)
    dc0 = K0 .* (1 .- exp.(-PRESSURE_B .* e)) ./ D0
    return c0, dc0
end

# Distance scaling of the registry-dependent stacking energy.
function registry_scale(d)
    e = d ./ D0 .- 1
    R = exp.(evalpoly_desc(LOGR, e))
    dR = R .* evalpoly_desc(derivpoly_desc(LOGR), e) ./ D0
    return R, dR
end

# Apply the elastic operator H (Fourier multiplier) to u; returns H u in real space.
function elastic_operator(h::VerticalHull, u)
    uk = fft(u, (2, 3, 4, 5))
    Hu = zeros(ComplexF64, size(uk))
    for a in 1:2, b in 1:2, c in 1:2, d in 1:2
        coef = h.C4[a, b, c, d]
        coef == 0 && continue
        @views Hu[a, :, :, :, :, :] .+= coef .* h.tK[b, :, :, :, :, :] .*
                                        h.tK[d, :, :, :, :, :] .* uk[c, :, :, :, :, :]
    end
    return real.(ifft(Hu, (2, 3, 4, 5)))
end

# Real-space Laplacian of each layer's field w (N,N,N,N,3).
function laplacian(h::VerticalHull, w)
    k2 = dropdims(sum(h.tK .^ 2; dims = 1); dims = 1)
    return real.(ifft(-k2 .* fft(w, (1, 2, 3, 4)), (1, 2, 3, 4)))
end

# Energy and gradient; x = [vec(u); vec(w); dbar12; dbar23].
function energy_gradient!(G, x, h::VerticalHull)
    N = h.N
    M = N^4
    nu = 2 * M * 3
    nw = M * 3
    u = reshape(view(x, 1:nu), (2, N, N, N, N, 3))
    w = reshape(view(x, nu+1:nu+nw), (N, N, N, N, 3))
    dbar = x[nu+nw+1:nu+nw+2]
    gu = zeros(2, N, N, N, N, 3)
    gw = zeros(N, N, N, N, 3)
    gd = zeros(2)

    Hu = elastic_operator(h, collect(u))
    F = 0.5 * dot(u, Hu) / M
    gu .+= Hu ./ M

    Lw = laplacian(h, collect(w))
    F += 0.5 * h.κ * h.Ω * sum(abs2, Lw) / M
    gw .+= h.κ * h.Ω .* laplacian(h, Lw) ./ M

    u2 = u[:, :, :, :, :, 2]
    w2 = w[:, :, :, :, 2]
    Fmis = 0.0
    Fbind = 0.0
    for (p, L, perm, sgn) in ((1, 1, h.perm1, -1.0), (2, 3, h.perm2, +1.0))
        uL = u[:, :, :, :, :, L]
        u2p = u2[:, perm]                                  # (2,N,N,N,N)
        dloc = dbar[p] .+ sgn .* (w[:, :, :, :, L] .- w2[perm])
        c0, dc0 = binding(h, dloc)
        R, dR = registry_scale(dloc)
        diff = u2p .- uL
        Gs = zeros(N, N, N, N)
        for invE in (h.invtE[2], h.invtE[L])
            shift = h.Γ0 .+ reshape(invE * reshape(diff, 2, M), (2, N, N, N, N))
            flat = reshape(shift, 2, M)
            Gs .+= 0.5 .* reshape(GSFE(flat), (N, N, N, N))
            gflat = gradient_GSFE(flat) .* (0.5 .* reshape(R, 1, M) ./ M)
            gdiff = reshape(invE' * gflat, (2, N, N, N, N))
            gu[:, :, :, :, :, L] .-= gdiff
            for I in CartesianIndices((N, N, N, N))
                gu[:, perm[I], 2] .+= gdiff[:, I]
            end
        end
        Fmis += sum(R .* Gs) / M
        Fbind += sum(c0) / M
        gdist = (dc0 .+ dR .* Gs) ./ M
        gw[:, :, :, :, L] .+= sgn .* gdist
        for I in CartesianIndices((N, N, N, N))
            gw[perm[I], 2] -= sgn * gdist[I]
        end
        gd[p] = sum(gdist) + h.pressure * h.Ω * GPaA3
    end
    F += Fmis + Fbind + h.pressure * h.Ω * GPaA3 * sum(dbar)
    for L in 1:3                                           # zero-mean gauge of w
        gw[:, :, :, :, L] .-= sum(gw[:, :, :, :, L]) / M
    end
    if G !== nothing
        G[1:nu] .= vec(gu)
        G[nu+1:nu+nw] .= vec(gw)
        G[nu+nw+1:end] .= gd
    end
    return F
end
