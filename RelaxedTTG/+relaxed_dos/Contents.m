% RELAXED_DOS  Paper DoS observables for relaxed, non-driven H(q).
%
% Public entry points
%   computeMomentumLocalDOS        - six-center-orbital momentum LDoS.
%   computeMomentumLocalDOSLinecut - momentum LDoS on a symmetry path.
%   computeTotalMomentumDOS        - reciprocal-space total DoS.
%   computeMoments                 - shared KPM moment/batching stage.
%   attachRelaxedHamiltonianToDoS  - attach cached relaxed H(q) to a Job.
%   plotMomentumLocalDOSLinecut    - plot the line-cut result.
%
% This module uses the paper's Jackson-Chebyshev regularization. It does
% not implement Green-operator integrated densities and contains no driven
% Hamiltonian routines.
