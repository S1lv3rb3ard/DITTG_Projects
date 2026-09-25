function [result,bridge] = computeRelaxedTotalIDOS( ...
    stack,DoF,shellsReference,relaxationFields,qPoints,qWeights,E,eta,options)
%COMPUTERELAXEDTOTALIDOS Total Green-operator IDoS for relaxed H(q).
%
% The default selector contains the six (G=0,j,alpha) center orbitals. The
% supplied qWeights must come from a positive reciprocal-space quadrature.
if nargin < 9
    options = struct;
end
defaults = struct( ...
    'hamiltonian',struct, ...
    'green',struct, ...
    'matrixCacheSize',1, ...
    'selector',[], ...
    'verbose',true);
options = mergeOptions(options,defaults);

if size(qPoints,1) ~= 2 && size(qPoints,2) == 2
    qPoints = qPoints.';
end
assert(size(qPoints,1) == 2 && ~isempty(qPoints), ...
    'qPoints must be a nonempty 2-by-N or N-by-2 array.');

cache = prepareRelaxedHamiltonian( ...
    stack,DoF,shellsReference,relaxationFields,qPoints,options.hamiltonian);
Hfun = rttg_common.makeCachedRelaxedHamiltonianFunction( ...
    cache,options.matrixCacheSize);
if isempty(options.selector)
    selector = getCenterBasis(DoF,1:6);
else
    selector = options.selector;
end

greenOptions = options.green;
greenOptions.verbose = options.verbose;
result = relaxed_idos.computeTotalIDOS( ...
    Hfun,qPoints,qWeights,selector,E,eta,greenOptions);

bridge.cache = cache;
bridge.H = Hfun;
bridge.selector = selector;
bridge.qPoints = qPoints;
bridge.qWeights = result.qWeights;
end
