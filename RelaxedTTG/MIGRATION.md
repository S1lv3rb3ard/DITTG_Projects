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
