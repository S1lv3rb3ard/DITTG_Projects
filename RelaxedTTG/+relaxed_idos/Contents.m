% RELAXED_IDOS  Integrated spectral measures for relaxed, non-driven H(q).
%
% Public entry points
%   computeMomentumLocalIDOS                 - retain q dependence.
%   computeRelaxedMomentumLocalIDOSLinecut   - symmetry-line wrapper.
%   computeTotalIDOS                         - positive q quadrature.
%   computeRelaxedTotalIDOS                  - relaxed-H(q) wrapper.
%   estimateFractalDimensionsFromIDOS        - finite-scale analysis.
%   plotMomentumLocalIDOSLinecut             - cumulative surface plot.
%
% This module uses Green operators, Poisson smoothing, and block Lanczos.
% It does not call the Jackson-KPM DoS module and contains no driven code.
