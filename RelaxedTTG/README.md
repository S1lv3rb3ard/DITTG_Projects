# Relaxed TTG Hamiltonian and spectral modules

This project assembles the relaxed, double-incommensurate trilayer
Hamiltonian

\[
H(q)=H_{\mathrm{intra}}(q)+H_{\mathrm{inter}}(q)
\]

and exposes three deliberately separate MATLAB modules. Run

```matlab
setupRelaxedTTG;
```

once before using the package-qualified APIs.

## Module boundaries

| MATLAB namespace | Hamiltonian | Observable and regularization |
|---|---|---|
| `relaxed_dos.*` | relaxed, non-driven `H(q)` | paper momentum LDoS and total DoS; Jackson-Chebyshev/KPM |
| `relaxed_idos.*` | relaxed, non-driven `H(q)` | momentum-local and reciprocal-space total IDoS; Green operator and Poisson kernel |
| `driven_rttg.*` | relaxed opto-acoustic space-time operator | central-sector driven spectral observables |

The three namespaces do not call one another. Shared numerical machinery is
in `rttg_common.*`; it defines no physical observable. Static relaxed
Hamiltonian assembly remains in `core/`, geometry in `geometry/`, hopping
models in `hopping/`, and general quadrature/device helpers in `utils/`.

## Relaxed Hamiltonian

The direct entry point is

```matlab
[H,parts,info] = buildRelaxedHamiltonian( ...
    stack,DoF,q,shells,relaxationFields,options);
```

`parts.intralayer` and `parts.interlayer` contain the two contributions.
For layer `j`, define `relaxationFields{j}(xk,xl)` with
`[k,l]=setdiff(1:3,j,'stable')`. Inputs and outputs are `2-by-N` physical
coordinate arrays: `xk` and `xl` are the position of a point of layer `j`
modulo the lattices of layers `k` and `l`, and the output is its in-plane
displacement in Angstrom.

### Relaxation fields

Mechanically relaxed fields come from the Julia configuration-space
minimizer (`Trilayers.jl`, `GrapheneParameters.jl`, `example2.jl`):

```bash
julia example2.jl <theta12> <theta23> <N>   # writes data/triG_*_<t12>_<t23>_<N>.jld
```

```matlab
stack = getStack(a,[-theta12,0,-theta23]);   % Julia rotates clockwise
[relaxationFields,info] = getRelaxationFields(stack,fullfile(projectRoot,'data'));
% or explicitly
[relaxationFields,info] = loadJuliaTrilayerRelaxation(minimizerFile,dataFile,stack);
```

`loadJuliaTrilayerRelaxation` reads the raw `Base.write` output (the `.jld`
extension is historical; the files are plain little-endian binaries) and
builds an exactly periodic trigonometric interpolant of the hull data. It
reconciles three differences between the two codes: Julia stores the shift
`b = R - x` (minus the MATLAB configuration), orders layer 3's configuration
pairs as (layer 2, layer 1), and uses a lattice basis rotated by 30 degrees
with clockwise twists. The frame rotation is found and checked against every
layer and orbital, and a stack whose angles do not match the Julia run is
rejected. `getRelaxationFields` picks the largest `N` available for the
stack's angles; the examples fall back to `makeToyTrilayerRelaxation` with a
warning when no Julia output is present.

**Evaluation points.** An orbital at `R + tau` in layer `j` is displaced by
`u_j(T_j R + D_j tau)`, where `D_j tau = (I - A_t A_j^{-1}) tau` in each slot
`t ~= j` (`getDisregistryOffsets`). This is the smooth disregistry of that
orbital. Evaluating at `T_j (R + tau)` instead shifts the configuration by a
sizeable fraction of a moire period and gives spurious bond strains of the
order of `max|u|` (for a (0.8, 1.1) degree relaxation, 0.27 A on
nearest-neighbour A-B bonds, against a physical 0.004 A). All intralayer
channels are sampled by `sampleRelaxedIntralayerChannels`, and interlayer
bonds by `evaluateRelaxedInterlayerPosition`.

**Vertical compression.** `getStack` sets both interfaces to `d0 = 3.35 A`.
At zero compression the interlayer hopping is exactly the Fang-Kaxiras
(2016) model used in arXiv:2606.13434 (model `'fangKaxiras2016Carr'`, the
default). Under compression or expansion, `eps = d/d0 - 1 ~= 0` (negative
under compression), the parameters follow the fits of Carr, Fang,
Jarillo-Herrero & Kaxiras (PRB 98, 085144, 2018) *relative to* `eps = 0`:
amplitudes are multiplied by `lambda_C(eps)/lambda_C(0)` and shape
parameters are shifted by `y_C(eps) - y_C(0)`.

```matlab
stack = setInterlayerCompression(stack,[eps12,eps23]);
stack = setInterlayerCompression(stack,compressionFromPressure(P_GPa));
```

The fits cover `-0.2 <= eps <= 0`. For expansion, `0 < eps <= 0.15` (for
example the AA regions of a corrugated moire), the amplitudes continue
exponentially and the shape parameters linearly, matched in value and slope
at `eps = 0`. Values outside `[-0.2, 0.15]` warn. A per-sample distance can
be passed directly,
`realSpaceInterlayerHopping(r,stack,j,k,alpha,beta,localEpsilon)`.
The raw Carr et al. fits are available as the model `'carr2018'`; note that
at `eps = 0` they differ from Fang-Kaxiras by 2-3% in the leading
tunnelling amplitude (and ~20% in `xi6`, `x6`), which visibly moves flat
bands near magic angles.

**Tests.** Run `runModuleTests("all")`. Local test functions must begin
with `test` to be collected by `functiontests`.
`tests/common/testUnrelaxedLimitMatchesPaper.m` checks that zero relaxation
reproduces the momentum-space Hamiltonian of arXiv:2606.13434.

### Corrugated, pressure-consistent relaxation (Julia)

`TrilayersVertical.jl` / `example_vertical.jl` extend the in-plane
minimizer with out-of-plane fields `w_j`, mean spacings `dbar12`, `dbar23`,
bending energy `(1/2) kappa Omega <(Lap w)^2>` (kappa = 1.4 eV by default),
the enthalpy `P Omega (dbar12 + dbar23)`, and a distance-dependent stacking
energy

    Phi(s,d) = c0(d) + R(d) * GSFE(s).

`GSFE` is the zero-pressure functional of `GrapheneParameters.jl`; `R(d)`
is the distance scaling of the AA-AB energy of the refined Kolmogorov-Crespi
potential (Ouyang et al., Nano Lett. 18, 6009 (2018)), fitted as
`log R = poly(eps)` on `-0.2 <= eps <= 0.15`; `c0(d)` integrates the pressure
law of Carr et al., so zero pressure gives `d0 = 3.35 A` for the
registry-averaged stack and the pressure scale matches the hopping fits.

```bash
julia example_vertical.jl <theta12> <theta23> <N> [P_GPa] [kappa_eV]
# -> data/triG_vz_{data,minimizer,vertical}_<t12>_<t23>_<N>_P<P>.jld
```

The `data` and `minimizer` files have the same layout as `example2.jl`
output, so `loadJuliaTrilayerRelaxation` reads the in-plane part unchanged.
Loading `w` into the interlayer hopping of RelaxedTTG is the next step.
The binding curve for expansion (extrapolated pressure law) and the
bending modulus are the main modelling uncertainties.

**Intralayer Bloch convention.** The intralayer blocks use the same Bloch
convention as the interlayer blocks: the element (G',j alpha),(G'',j beta)
is `sum_d [h]_{G''-G'}(d) exp(-i(q+SG').bond0) exp(i S(G''-G').tau_beta)` with
`bond0 = d + tau_alpha - tau_beta`, stored per channel in the cache as
`bondAA`, `bondAB`, `bondBA`, `bondBB`.

The interlayer transform uses no FFT. It evaluates the continuous
real-space Fourier integral by direct polar quadrature and the normalized
spectator-cell integral by a direct periodic trapezoidal sum. State pairs
with the same spectator transfer are grouped so each configuration
coefficient is integrated only once.

`options.interlayer.spectatorModeCutoff` retains transfers satisfying

```matlab
max(abs(n_l)) <= spectatorModeCutoff
```

where `G_l''-G_l'=B_l*n_l`. Use `inf` for the exact finite basis and converge
any finite value. Converge `W`, `L`, configuration grids, spectator cutoff,
real-space cutoff, and radial/angular quadrature independently.

Hopping-function visualization remains in `plotting/`. In particular,
`plotInterlayerMassattComparison` can use `transformType='spectatorSlice'`
to transform exactly the displayed fixed-spectator slice. Its Fourier
transform is continuous direct quadrature, not an FFT.

## Module 1: `relaxed_dos`

This module implements the paper's relaxed, non-driven density-of-states
observables. The momentum LDoS is exactly the six-center-orbital quantity

\[
\widehat{\mathcal D}_{\epsilon}(E,q)
=\frac{1}{6}\sum_{j\alpha}
 [\delta_{\epsilon}(E-\widehat H(q))]_{(0,j\alpha),(0,j\alpha)}.
\]

It preserves the established selector and normalization rather than
redefining the observable through the Green-IDoS code.

Public entry points are

```matlab
Job = relaxed_dos.computeMomentumLocalDOS(Job,E,mode);
Job = relaxed_dos.computeTotalMomentumDOS(Job,E,mode);

[Job,bridge] = relaxed_dos.attachRelaxedHamiltonianToDoS( ...
    Job,thetaIndex,stack,DoF,shells,relaxationFields,options);

Job = relaxed_dos.computeMoments(Job,mode);
Job = relaxed_dos.computeMomentumLocalDOS(Job,E,mode); % or total

[result,Job,bridge] = relaxed_dos.computeMomentumLocalDOSLinecut( ...
    stack,DoF,shells,relaxationFields,linecut,E,options);
```

Construct and plot the paper's line cut with

```matlab
linecut = rttg_common.getLinecut(stack.K,pointsPerSegment);
relaxed_dos.plotMomentumLocalDOSLinecut( ...
    linecut,E,result.LDoS);
```

The high-symmetry example is

```matlab
run('examples/relaxed_dos/runMomentumLocalDOSLinecut.m')
```

`P` is the Chebyshev order. A value such as `P=256` is generally only a
code-path check for the energy scales of interest; the included diagnostic
example uses `P=4000` and reports the central Jackson-kernel resolution.

## Module 2: `relaxed_idos`

This module implements cumulative spectral measures through

\[
G_q(z)=(H(q)-zI)^{-1},\qquad z=E+i\eta,
\]

using block Lanczos with two-pass orthogonal correction. Its two distinct
public calculations are

```matlab
local = relaxed_idos.computeMomentumLocalIDOS( ...
    Hamiltonian,qPoints,selector,E,eta,options);

total = relaxed_idos.computeTotalIDOS( ...
    Hamiltonian,qPoints,qWeights,selector,E,eta,options);
```

The local call retains one result for every q-point:

```matlab
local.momentumLocalIDOS(qIndex,energyIndex,etaIndex)
```

The total call applies a positive normalized reciprocal-cell quadrature and
returns

```matlab
total.totalIDOS(etaIndex,energyIndex)
```

It has no remaining q-axis. The corresponding relaxed-Hamiltonian wrappers
are

```matlab
[local,bridge] = ...
    relaxed_idos.computeRelaxedMomentumLocalIDOSLinecut( ...
        stack,DoF,shells,relaxationFields,linecut,E,eta,options);

[total,bridge] = relaxed_idos.computeRelaxedTotalIDOS( ...
    stack,DoF,shells,relaxationFields, ...
    qPoints,qWeights,E,eta,options);
```

For the cumulative Poisson regularization,

\[
N_\eta(q,E)=\int
\left[\frac12+\frac1\pi\arctan\frac{E-\lambda}{\eta}\right]
\,d\mu_q(\lambda).
\]

The paper's Jackson-KPM momentum LDoS is not returned from this module; it
belongs exclusively to `relaxed_dos`.

Plot a cumulative high-symmetry surface with

```matlab
relaxed_idos.plotMomentumLocalIDOSLinecut( ...
    linecut,E,local.momentumLocalIDOS,eta);
```

Examples are

```matlab
run('examples/relaxed_idos/runMomentumLocalIDOSLinecut.m')
run('examples/relaxed_idos/runTotalIDOSFractalDimension.m')
```

The total example uses a positive reciprocal-cell trapezoidal quadrature;
it is deliberately coarse and must be converged. Fractal analysis uses

```matlab
dimension = relaxed_idos.estimateFractalDimensionsFromIDOS( ...
    E,total.totalIDOS,incrementRadii,options);
```

Do not send `eta` to zero at fixed finite truncation or fixed Lanczos depth.
That limit only resolves the artificial atomic spectrum of the finite
matrix. Converge spatial truncation, reciprocal quadrature, Lanczos depth,
`eta`, and the IDoS increment scale jointly.

## Module 3: `driven_rttg`

This module owns the doubly incommensurate opto-acoustic space-time
Hamiltonian and its driven observables. For optical sideband indices `n`
and acoustic indices `m`, diagonal blocks use

```text
H(q - n*K - m*Q) + hbar*(n*omega + m*Omega)*I.
```

Multiple optical and acoustic sources are supported. Optical sources use
repeating pulse trains and are checked against the configured duty-cycle
guard.

```matlab
optical(1) = driven_rttg.makeOpticalDriveSource( ...
    amplitude,K,omega,phase,cutoff,envelope);
acoustic(1) = driven_rttg.makeAcousticDriveSource( ...
    amplitude,Q,Omega,phase,cutoff);

[driveCache,bridge] = ...
    driven_rttg.prepareRelaxedDrivenSpaceTimeHamiltonian( ...
        stack,DoF,shells,relaxationFields,q, ...
        optical,acoustic,options);

[K,parts,info] = driven_rttg.evaluateDrivenSpaceTimeHamiltonian( ...
    driveCache,envelopeTime);

observable = driven_rttg.computeDrivenSpectralObservables( ...
    driveCache,E,eta,envelopeTime,observableOptions);
```

The driven observable is a central-temporal-sector, frozen-envelope
spectral diagnostic. It is not the non-driven momentum LDoS or IDoS, and it
is not a time-resolved ARPES intensity. The latter requires propagation or
a lesser Green function together with occupations, probe envelopes, and
photoemission matrix elements.

Run

```matlab
run('examples/driven_rttg/runIncommensuratelyDrivenTTG.m')
```

for the multiple-source repeating-pulse example.

## Tests

Tests follow the same module boundaries:

```matlab
results = runModuleTests("relaxed_dos");
assertSuccess(results);

results = runModuleTests("relaxed_idos");
assertSuccess(results);

results = runModuleTests("driven_rttg");
assertSuccess(results);
```

Use `runModuleTests("all")` for the complete suite. The `common` tests cover
Hamiltonian and hopping conventions, not a fourth observable module.
