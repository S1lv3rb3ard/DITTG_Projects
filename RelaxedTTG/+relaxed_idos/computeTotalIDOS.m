function result = computeTotalIDOS( ...
    Hamiltonian,qPoints,qWeights,selector,E,eta,options)
%COMPUTETOTALIDOS Evaluate reciprocal-space total IDoS from Green operators.
%
% result = relaxed_idos.computeTotalIDOS( ...
%     Hamiltonian,qPoints,qWeights,selector,E,eta,options)
% runs fully reorthogonalized block Lanczos once at every q. Hamiltonian may
% be a matrix, a cell array of matrices, or a function H(q). selector may
% likewise be a matrix, a cell array, or a function V(q). qWeights must be
% nonnegative and are normalized internally.
%
% The convention is G(z)=(H-z I)^(-1), so
%
%   rho_eta(E) = Im m(E+i eta)/pi,
%   N_eta(E)   = integral_{-inf}^E rho_eta(lambda) dlambda.
%
% N_eta is evaluated analytically from the positive Lanczos-Pade weights,
% rather than by numerically integrating a sampled density.
if nargin < 7
    options = struct;
end
defaults = struct( ...
    'lanczos',struct, ...
    'energyChunkSize',512, ...
    'storeLanczos',false, ...
    'verbose',true);
options = mergeOptions(options,defaults);

[qPoints,numberQ] = normalizeProviders(Hamiltonian,qPoints);
assert(isvector(qWeights) && numel(qWeights) == numberQ, ...
    'qWeights must contain one entry per Hamiltonian.');
qWeights = qWeights(:);
assert(all(isfinite(qWeights)) && all(qWeights >= 0) && sum(qWeights) > 0, ...
    'qWeights must be finite, nonnegative, and have positive sum.');
qWeights = qWeights/sum(qWeights);
assert(isvector(E) && ~isempty(E) && all(isfinite(E)), ...
    'E must be a finite, nonempty vector.');
assert(isvector(eta) && ~isempty(eta) && ...
    all(isfinite(eta)) && all(eta > 0), ...
    'eta must be a finite vector of positive smoothing widths.');
assert(options.energyChunkSize >= 1 && ...
    options.energyChunkSize == fix(options.energyChunkSize), ...
    'energyChunkSize must be a positive integer.');

E = E(:).';
eta = eta(:);
allNodes = cell(numberQ,1);
allWeights = cell(numberQ,1);
diagnostics = repmat(struct,numberQ,1);
if options.storeLanczos
    recursion = cell(numberQ,1);
else
    recursion = {};
end

for qIndex = 1:numberQ
    Hq = evaluateProvider(Hamiltonian,qPoints,qIndex);
    Vq = evaluateProvider(selector,qPoints,qIndex);
    lanczos = rttg_common.blockLanczosHermitian( ...
        Hq,Vq,options.lanczos);
    allNodes{qIndex} = lanczos.ritzValues;
    allWeights{qIndex} = qWeights(qIndex)*lanczos.scalarWeights;
    diagnostics(qIndex).numberBlockSteps = lanczos.numberBlockSteps;
    diagnostics(qIndex).krylovDimension = lanczos.krylovDimension;
    diagnostics(qIndex).orthogonalityError = lanczos.orthogonalityError;
    diagnostics(qIndex).blockTridiagonalDefect = ...
        lanczos.blockTridiagonalDefect;
    diagnostics(qIndex).weightError = abs(sum(lanczos.scalarWeights)-1);
    if options.storeLanczos
        recursion{qIndex} = lanczos;
    end
    if options.verbose
        fprintf('Green IDoS: completed q-point %d of %d.\n',qIndex,numberQ);
    end
end

nodes = vertcat(allNodes{:});
weights = vertcat(allWeights{:});
weights = max(real(weights),0);
weights = weights/sum(weights);

numberEta = numel(eta);
numberEnergy = numel(E);
green = complex(zeros(numberEta,numberEnergy));
IDOS = zeros(numberEta,numberEnergy);
for firstEnergy = 1:options.energyChunkSize:numberEnergy
    lastEnergy = min(firstEnergy+options.energyChunkSize-1,numberEnergy);
    range = firstEnergy:lastEnergy;
    energyBlock = E(range);
    for etaIndex = 1:numberEta
        z = energyBlock+1i*eta(etaIndex);
        green(etaIndex,range) = weights.'*(1./(nodes-z));
        IDOS(etaIndex,range) = weights.'*( ...
            0.5+atan((energyBlock-nodes)/eta(etaIndex))/pi);
    end
end
assert(min(imag(green),[],'all') >= -1e-12, ...
    'The computed Green function violates the Herglotz sign condition.');

[sortedNodes,order] = sort(nodes);
sortedWeights = weights(order);
cumulativeWeights = cumsum(sortedWeights);
distinctNodes = unique(sortedNodes);
if numel(distinctNodes) >= 2
    poleSpacing = diff(distinctNodes);
    effectivePoleSpacing = median(poleSpacing);
    minimumPoleSpacing = min(poleSpacing);
else
    effectivePoleSpacing = inf;
    minimumPoleSpacing = inf;
end
padeIDOS = zeros(1,numberEnergy);
for energyIndex = 1:numberEnergy
    location = find(sortedNodes <= E(energyIndex),1,'last');
    if ~isempty(location)
        padeIDOS(energyIndex) = cumulativeWeights(location);
    end
end

result.energy = E;
result.eta = eta;
result.green = green;
result.totalIDOS = IDOS;
result.padeIDOS = padeIDOS;
result.padeNodes = nodes;
result.padeWeights = weights;
result.effectivePoleSpacing = effectivePoleSpacing;
result.minimumPoleSpacing = minimumPoleSpacing;
result.qPoints = qPoints;
result.qWeights = qWeights;
result.diagnostics = diagnostics;
result.lanczos = recursion;
result.convention = 'G(z)=(H-zI)^(-1), Im G(z)>0 for Im z>0';
result.warning = ['Finite Lanczos-Pade measures are atomic. Estimate ' ...
    'dimensions only in a joint size/depth/smoothing limit above the ' ...
    'effective pole spacing.'];
result.observable = 'reciprocal-space total integrated density of states';
end

function [qPoints,numberQ] = normalizeProviders(Hamiltonian,qPoints)
if isnumeric(Hamiltonian)
    numberQ = 1;
elseif iscell(Hamiltonian)
    numberQ = numel(Hamiltonian);
elseif isa(Hamiltonian,'function_handle')
    assert(~isempty(qPoints), ...
        'qPoints are required when Hamiltonian is a function handle.');
    if size(qPoints,1) ~= 2 && size(qPoints,2) == 2
        qPoints = qPoints.';
    end
    assert(size(qPoints,1) == 2,'qPoints must be 2-by-N or N-by-2.');
    numberQ = size(qPoints,2);
else
    error('Unsupported Hamiltonian provider.');
end
if isempty(qPoints)
    qPoints = zeros(2,numberQ);
elseif size(qPoints,1) ~= 2 && size(qPoints,2) == 2
    qPoints = qPoints.';
end
assert(size(qPoints,1) == 2 && size(qPoints,2) == numberQ, ...
    'qPoints must contain one column per Hamiltonian.');
end

function value = evaluateProvider(provider,qPoints,index)
if isnumeric(provider)
    value = provider;
elseif iscell(provider)
    assert(numel(provider) >= index, ...
        'The provider cell array has too few entries.');
    value = provider{index};
elseif isa(provider,'function_handle')
    if nargin(provider) == 1
        value = provider(qPoints(:,index));
    else
        value = provider(qPoints(:,index),index);
    end
else
    error('Provider must be numeric, a cell array, or a function handle.');
end
end
