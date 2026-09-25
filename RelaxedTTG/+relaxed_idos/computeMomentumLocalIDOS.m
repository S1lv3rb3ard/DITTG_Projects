function result = computeMomentumLocalIDOS( ...
    Hamiltonian,qPoints,selector,E,eta,options)
%COMPUTEMOMENTUMLOCALIDOS Momentum-local Green IDoS at retained q-points.
%
% result = relaxed_idos.computeMomentumLocalIDOS( ...
%     Hamiltonian,qPoints,selector,E,eta,options)
% preserves the q dependence instead of applying reciprocal-space
% quadrature.  Its output is the cumulative local spectral measure
%
%   N_eta(q,E) = sum_l w_l(q) [1/2 + atan((E-lambda_l(q))/eta)/pi].
%
% result.momentumLocalIDOS has size
% numberQ-by-numberEnergy-by-numberEta. For a single eta,
% squeeze(result.momentumLocalIDOS(:,:,1)) is the high-symmetry surface.
% This is a momentum-local integrated spectral measure, not the total IDoS
% obtained after integrating over q.
if nargin < 6
    options = struct;
end
defaults = struct( ...
    'lanczos',struct, ...
    'energyChunkSize',512, ...
    'storeLanczos',false, ...
    'verbose',true);
options = mergeOptions(options,defaults);

if size(qPoints,1) ~= 2 && size(qPoints,2) == 2
    qPoints = qPoints.';
end
assert(size(qPoints,1) == 2 && ~isempty(qPoints), ...
    'qPoints must be a nonempty 2-by-N or N-by-2 array.');
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
numberQ = size(qPoints,2);
numberEnergy = numel(E);
numberEta = numel(eta);

green = complex(zeros(numberQ,numberEnergy,numberEta));
IDOS = zeros(numberQ,numberEnergy,numberEta);
padeIDOS = zeros(numberQ,numberEnergy);
padeNodes = cell(numberQ,1);
padeWeights = cell(numberQ,1);
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

    nodes = real(lanczos.ritzValues(:));
    weights = max(real(lanczos.scalarWeights(:)),0);
    weights = weights/sum(weights);
    padeNodes{qIndex} = nodes;
    padeWeights{qIndex} = weights;

    diagnostics(qIndex).numberBlockSteps = lanczos.numberBlockSteps;
    diagnostics(qIndex).krylovDimension = lanczos.krylovDimension;
    diagnostics(qIndex).orthogonalityError = lanczos.orthogonalityError;
    diagnostics(qIndex).blockTridiagonalDefect = ...
        lanczos.blockTridiagonalDefect;
    diagnostics(qIndex).weightError = abs(sum(weights)-1);
    if options.storeLanczos
        recursion{qIndex} = lanczos;
    end

    for firstEnergy = 1:options.energyChunkSize:numberEnergy
        lastEnergy = min(firstEnergy+options.energyChunkSize-1,numberEnergy);
        range = firstEnergy:lastEnergy;
        energyBlock = E(range);
        for etaIndex = 1:numberEta
            z = energyBlock+1i*eta(etaIndex);
            green(qIndex,range,etaIndex) = ...
                weights.'*(1./(nodes-z));
            IDOS(qIndex,range,etaIndex) = weights.'*( ...
                0.5+atan((energyBlock-nodes)/eta(etaIndex))/pi);
        end
    end

    [sortedNodes,order] = sort(nodes);
    cumulativeWeights = cumsum(weights(order));
    for energyIndex = 1:numberEnergy
        location = find(sortedNodes <= E(energyIndex),1,'last');
        if ~isempty(location)
            padeIDOS(qIndex,energyIndex) = cumulativeWeights(location);
        end
    end

    if options.verbose
        fprintf('Green line cut: completed q-point %d of %d.\n', ...
            qIndex,numberQ);
    end
end

assert(min(imag(green),[],'all') >= -1e-12, ...
    'The computed Green function violates the Herglotz sign condition.');

result.energy = E;
result.eta = eta;
result.qPoints = qPoints;
result.green = green;
result.momentumLocalIDOS = IDOS;
result.padeIDOS = padeIDOS;
result.padeNodes = padeNodes;
result.padeWeights = padeWeights;
result.diagnostics = diagnostics;
result.lanczos = recursion;
result.convention = 'G(z)=(H-zI)^(-1), Im G(z)>0 for Im z>0';
result.interpretation = [ ...
    'Momentum-local cumulative spectral measure. Integrate over a ' ...
    'positive reciprocal-space quadrature to obtain the total IDoS.'];
result.observable = 'momentum-local integrated density of states';
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
