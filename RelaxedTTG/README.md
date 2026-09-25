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
coordinate arrays.

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
