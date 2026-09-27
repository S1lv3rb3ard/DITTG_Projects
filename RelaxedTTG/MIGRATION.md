# Module migration map

The previous flat folders mixed non-driven DoS, non-driven IDoS, shared
Lanczos machinery, and driven observables. Use the following replacements.

| Previous call | New call |
|---|---|
| `getLDoS(...)` | `relaxed_dos.computeMomentumLocalDOS(...)` |
| `getDoS(...)` | `relaxed_dos.computeTotalMomentumDOS(...)` |
| `attachRelaxedHamiltonianToDoS(...)` | `relaxed_dos.attachRelaxedHamiltonianToDoS(...)` |
| `BatchWeights(...)` | `relaxed_dos.computeMoments(...)` |
| `computeRelaxedMomentumLDoSLinecut(...)` | `relaxed_dos.computeMomentumLocalDOSLinecut(...)` |
| `plotMomentumLDoSLinecut(...)` | `relaxed_dos.plotMomentumLocalDOSLinecut(...)` |
| `computeGreenIDOSLinecut(...)` | `relaxed_idos.computeMomentumLocalIDOS(...)` |
| `computeRelaxedMomentumIDOSLinecut(...)` | `relaxed_idos.computeRelaxedMomentumLocalIDOSLinecut(...)` |
| `computeGreenIDOS(...)` | `relaxed_idos.computeTotalIDOS(...)` |
| `computeRelaxedGreenIDOS(...)` | `relaxed_idos.computeRelaxedTotalIDOS(...)` |
| `plotMomentumIDOSLinecut(...)` | `relaxed_idos.plotMomentumLocalIDOSLinecut(...)` |
| `estimateFractalDimensionsFromIDOS(...)` | `relaxed_idos.estimateFractalDimensionsFromIDOS(...)` |
| any former `driving/` call | prefix the function with `driven_rttg.` |

The IDoS result fields are now explicit:

```matlab
local.momentumLocalIDOS
total.totalIDOS
```

The former public `poissonDOS` field was removed from the IDoS API. The
paper's regularized momentum density of states is computed only by
`relaxed_dos`; this prevents a Poisson derivative of an IDoS approximation
from being silently substituted for the paper's Jackson-KPM observable.

Shared helpers moved to `rttg_common.*`:

```matlab
rttg_common.blockLanczosHermitian
rttg_common.evaluateBlockLanczosGreen
rttg_common.getLinecut
rttg_common.makeCachedRelaxedHamiltonianFunction
```

## Relaxation evaluation points and intralayer gauge (2026-09)

* Orbital displacements are now evaluated at the disregistry configuration
  `T_j R + D_j tau` (`getDisregistryOffsets`) instead of `T_j (R + tau)`. The
  intralayer AB/BA and BB channels and all interlayer channels with a B
  orbital change; AA channels are unchanged.
* The intralayer A-B blocks previously mixed two Bloch gauges (the diagonal
  symmetrization averaged `sum t exp(-iq.(RB-tau))` with the conjugate of
  `sum t exp(-iq.(RB+tau))`), which made the unrelaxed monolayer Dirac cone
  anisotropic, with the velocity halved along some directions. They now use
  the interlayer convention; cache fields `RA`, `RBminusTau`, `RBplusTau` are
  replaced by `bondAA`, `bondAB`, `bondBA`, `bondBB`.
* New: `loadJuliaTrilayerRelaxation`, `readJuliaRelaxationFiles`,
  `evaluateFourierRelaxationField`, `getRelaxationFields`,
  `sampleRelaxedIntralayerChannels`, `getDisregistryOffsets`, and the test
  `tests/common/testRelaxationConventions.m` with fixture `tests/data/`.

## Uniform interlayer compression (2026-09)

* `realSpaceInterlayerHopping` now takes its ten parameters from
  `stack.interlayer` (set by `getStack` and `setInterlayerCompression`). The
  default is the Carr et al. (2018) fit at zero compression, which differs
  slightly from the constants used before (for example xi3 = 3.286 instead of
  3.4692). Use `setInterlayerCompression(stack,0,'fangKaxiras2016')` to
  reproduce earlier spectra exactly.
* Only the adjacent interfaces (1,2) and (2,3) are parameterized.

## Expansion continuation, local distances, corrugated relaxation (2026-09)

* `getInterlayerHoppingParameters('carr2018',eps)` accepts arrays and, for
  `0 < eps <= 0.15`, uses a C^1 continuation (exponential amplitudes, linear
  shape parameters) instead of extrapolating the quadratic fits.
* `realSpaceInterlayerHopping` accepts an optional per-sample `localEpsilon`.
* New Julia files `TrilayersVertical.jl`, `example_vertical.jl` for
  corrugated, pressure-consistent relaxation.
