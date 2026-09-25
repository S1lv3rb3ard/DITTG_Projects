function [result,bridge] = computeRelaxedMomentumLocalIDOSLinecut( ...
    stack,DoF,shellsReference,relaxationFields,linecut,E,eta,options)
%COMPUTERELAXEDMOMENTUMLOCALIDOSLINECUT Local IDoS over a symmetry path.
%
% The result is momentum local:
%
%   result.momentumLocalIDOS(qIndex,energyIndex,etaIndex) = N_eta(q,E).
%
% It is not the total IDoS, which would require a positive quadrature over
% the reciprocal cell and would no longer retain a high-symmetry-path axis.
if nargin < 8
    options = struct;
end
defaults = struct( ...
    'hamiltonian',struct, ...
    'green',struct, ...
    'matrixCacheSize',1, ...
    'selector',[], ...
    'verbose',true);
options = mergeOptions(options,defaults);

assert(isstruct(linecut) && all(isfield(linecut, ...
    {'q','kPath','tickLocs','labels'})), ...
    'linecut must be produced by getLinecut.');
assert(size(linecut.q,1) == 2 && ...
    size(linecut.q,2) == numel(linecut.kPath), ...
    'linecut.q must be 2-by-N and agree with linecut.kPath.');

cache = prepareRelaxedHamiltonian(stack,DoF,shellsReference, ...
    relaxationFields,linecut.q,options.hamiltonian);
Hfun = rttg_common.makeCachedRelaxedHamiltonianFunction( ...
    cache,options.matrixCacheSize);
if isempty(options.selector)
    selector = getCenterBasis(DoF,1:6);
else
    selector = options.selector;
end

greenOptions = options.green;
greenOptions.verbose = options.verbose;
result = relaxed_idos.computeMomentumLocalIDOS( ...
    Hfun,linecut.q,selector,E,eta,greenOptions);
result.q = linecut.q;
result.kPath = linecut.kPath;
result.tickLocs = linecut.tickLocs;
result.labels = linecut.labels;
result.projector = selector;

bridge.cache = cache;
bridge.H = Hfun;
bridge.selector = selector;
bridge.qPoints = linecut.q;
end
