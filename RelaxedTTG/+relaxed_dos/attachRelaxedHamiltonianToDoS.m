function [Job,bridge] = attachRelaxedHamiltonianToDoS( ...
    Job,thetaIndex,stack,DoF,shellsReference,relaxationFields,options)
%ATTACHRELAXEDHAMILTONIANTODOS Connect cached relaxed H(q) to a DoS Job.
%
% This function preserves the existing momentum-LDoS and total-DoS
% projectors, KPM recursion, quadrature weights, and normalization.  It only
% replaces Job.H{thetaIndex} with the relaxed Hamiltonian constructor.
if nargin < 7
    options = struct;
end
defaults = struct( ...
    'hamiltonian',struct, ...
    'matrixCacheSize',1, ...
    'validateSelectors',true, ...
    'setMomentumLDoSSelector',false, ...
    'momentumLDoSSelectorIndex',1, ...
    'verbose',true);
options = mergeOptions(options,defaults);

assert(istable(Job.list),'Job.list must be a table.');
assert(all(ismember({'t','Q'},Job.list.Properties.VariableNames)), ...
    'Job.list must contain variables t and Q.');
assert(thetaIndex >= 1 && thetaIndex == fix(thetaIndex), ...
    'thetaIndex must be a positive integer.');
rows = Job.list.t == thetaIndex;
assert(any(rows),'Job.list contains no rows with t=%d.',thetaIndex);
qPoints = unique(Job.list.Q(rows,:),'rows','stable').';

cache = prepareRelaxedHamiltonian( ...
    stack,DoF,shellsReference,relaxationFields,qPoints,options.hamiltonian);
Hfun = rttg_common.makeCachedRelaxedHamiltonianFunction( ...
    cache,options.matrixCacheSize);

if ~isfield(Job,'H') || isempty(Job.H)
    Job.H = cell(1,thetaIndex);
elseif numel(Job.H) < thetaIndex
    Job.H{thetaIndex} = [];
end
Job.H{thetaIndex} = Hfun;

if options.setMomentumLDoSSelector
    if ~isfield(Job,'X') || isempty(Job.X)
        Job.X = cell(1,thetaIndex);
    elseif numel(Job.X) < thetaIndex
        Job.X{thetaIndex} = {};
    end
    if isempty(Job.X{thetaIndex})
        Job.X{thetaIndex} = cell(1,options.momentumLDoSSelectorIndex);
    elseif numel(Job.X{thetaIndex}) < options.momentumLDoSSelectorIndex
        Job.X{thetaIndex}{options.momentumLDoSSelectorIndex} = [];
    end
    Job.X{thetaIndex}{options.momentumLDoSSelectorIndex} = ...
        getCenterBasis(DoF,1:6);
end

if options.validateSelectors && isfield(Job,'X') && ...
        numel(Job.X) >= thetaIndex && ~isempty(Job.X{thetaIndex})
    selectors = Job.X{thetaIndex};
    for index = 1:numel(selectors)
        if isempty(selectors{index})
            continue
        end
        assert(size(selectors{index},1) == cache.dimension, ...
            ['Job.X{%d}{%d} has %d rows, but the relaxed Hamiltonian ' ...
             'dimension is %d. Rebuild the selector from the same DoF.'], ...
            thetaIndex,index,size(selectors{index},1),cache.dimension);
    end
end

bridge.cache = cache;
bridge.H = Hfun;
bridge.qPoints = qPoints;
bridge.thetaIndex = thetaIndex;
bridge.dimension = cache.dimension;
bridge.momentumLDoSBasis = getCenterBasis(DoF,1:6);
bridge.note = ['Job.scale and Job.weights are intentionally preserved. ' ...
    'Verify scale*H has spectrum in (-1,1), and build weights from the ' ...
    'same DoF truncation.'];

if options.verbose
    fprintf(['Attached relaxed H(q) to Job.H{%d}: %d q-points, ' ...
        'dimension %d.\n'],thetaIndex,size(qPoints,2),cache.dimension);
end
end
