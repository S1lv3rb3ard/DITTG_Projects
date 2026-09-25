function cache = prepareRelaxedHamiltonian( ...
    stack,DoF,shellsReference,relaxationFields,qPoints,options)
%PREPARERELAXEDHAMILTONIAN Cache q-independent relaxed Hamiltonian data.
%
% cache = prepareRelaxedHamiltonian(...,qPoints,options) prepares a reusable
% Hamiltonian cache for every q inside the smallest enclosing ball used by
% qPoints.  The expensive relaxation-field sampling, continuous real-space
% quadrature, and spectator-cell trapezoidal integration are performed here.
% evaluateRelaxedHamiltonian(cache,q) then performs only q-dependent phase,
% cutoff, and sparse-matrix assembly work.
%
% qPoints may be 2-by-N or N-by-2.  Supplying the actual DoS/LDoS sampling
% points gives the tightest interlayer candidate list.
if nargin < 6
    options = struct;
end
defaults = struct( ...
    'intralayer',struct, ...
    'interlayer',struct, ...
    'qPadding',0, ...
    'hermiticityTolerance',1e-12, ...
    'verbose',true);
options = mergeOptions(options,defaults);

qPoints = normalizeQPoints(qPoints);
assert(~isempty(qPoints),'qPoints must contain at least one point.');
assert(options.qPadding >= 0,'options.qPadding must be nonnegative.');
qReference = mean(qPoints,2);
qRadius = max(vecnorm(qPoints-qReference,2,1))+options.qPadding;

cache.stack = stack;
cache.DoF = DoF;
cache.dimension = 2*size(DoF,1);
cache.qReference = qReference;
cache.qRadius = qRadius;
cache.qPoints = qPoints;
cache.options = options;

cache.intralayer = prepareRelaxedIntralayerHamiltonian( ...
    stack,DoF,shellsReference,relaxationFields,options.intralayer);
cache.interlayer = prepareRelaxedInterlayerHamiltonian( ...
    stack,DoF,relaxationFields,qReference,qRadius,options.interlayer);

if options.verbose
    fprintf(['Prepared relaxed Hamiltonian cache for %d q-points ' ...
        '(domain radius %.6g, dimension %d).\n'], ...
        size(qPoints,2),qRadius,cache.dimension);
end
end

function qPoints = normalizeQPoints(qPoints)
assert(isnumeric(qPoints) && ismatrix(qPoints), ...
    'qPoints must be a numeric 2-by-N or N-by-2 array.');
if size(qPoints,1) == 2
    return
end
if size(qPoints,2) == 2
    qPoints = qPoints.';
    return
end
error('qPoints must be a numeric 2-by-N or N-by-2 array.');
end
