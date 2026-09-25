function cache = prepareRelaxedInterlayerHamiltonian( ...
    stack,DoF,relaxationFields,qReference,qRadius,options)
%PREPARERELAXEDINTERLAYERHAMILTONIAN Cache direct-quadrature interlayer data.
%
% No FFT is used.  The continuous R^2 samples and normalized spectator-cell
% trapezoidal coefficients are calculated once and retained for reuse.
if nargin < 6
    options = struct;
end
defaults = struct( ...
    'realSpaceRadialOrder',80, ...
    'realSpaceAngularOrder',128, ...
    'realSpaceCutoff',8*stack.aG, ...
    'configurationGridSize',[], ...
    'minimumConfigurationGridSize',9, ...
    'spectatorModeCutoff',inf, ...
    'momentumInnerRadius',6.00, ...
    'momentumOuterRadius',6.25, ...
    'modeBatchSize',16, ...
    'pairBatchSize',128, ...
    'useGPU',false, ...
    'verbose',true);
options = mergePreservingEmpty(options,defaults,'configurationGridSize');
validateOptions(options);

assert(iscell(relaxationFields) && numel(relaxationFields) == 3, ...
    'relaxationFields must contain three function handles.');
assert(isequal(size(qReference),[2,1]), ...
    'qReference must be a 2-by-1 column vector.');
assert(isscalar(qRadius) && qRadius >= 0,'qRadius must be nonnegative.');

cache.dimension = 2*size(DoF,1);
cache.qReference = qReference;
cache.qRadius = qRadius;
cache.options = options;
cache.pair = cell(1,2);
pairLayers = [1,2,3;2,3,1];
for pair = 1:2
    cache.pair{pair} = preparePair(stack,DoF,relaxationFields, ...
        qReference,qRadius,options,pairLayers(pair,1), ...
        pairLayers(pair,2),pairLayers(pair,3));
end
end

function data = preparePair( ...
    stack,DoF,fields,qReference,qRadius,options,j,k,l)
indexJ = find(DoF(:,9) == j);
indexK = find(DoF(:,9) == k);
assert(~isempty(indexJ) && ~isempty(indexK), ...
    'Layers %d and %d must contain retained DoFs.',j,k);

GkPrime = DoF(indexJ,2*k+(1:2));
GlPrime = DoF(indexJ,2*l+(1:2));
GjDouble = DoF(indexK,2*j+(1:2));
GlDouble = DoF(indexK,2*l+(1:2));
nJ = numel(indexJ);
nK = numel(indexK);
candidateRadius = options.momentumOuterRadius+qRadius;

pairCapacity = 0;
for sourceIndex = 1:nJ
    Q = qReference.'+GjDouble+GkPrime(sourceIndex,:)+GlPrime(sourceIndex,:);
    pairCapacity = pairCapacity+nnz(vecnorm(Q,2,2) < candidateRadius);
end

source = zeros(pairCapacity,1);
target = zeros(pairCapacity,1);
baseQ = zeros(pairCapacity,2);
dGl = zeros(pairCapacity,2);
ptr = 0;
for sourceIndex = 1:nJ
    Q = qReference.'+GjDouble+GkPrime(sourceIndex,:)+GlPrime(sourceIndex,:);
    retainedTarget = find(vecnorm(Q,2,2) < candidateRadius);
    number = numel(retainedTarget);
    if number == 0
        continue
    end
    range = ptr+(1:number);
    source(range) = sourceIndex;
    target(range) = retainedTarget;
    baseQ(range,:) = GjDouble(retainedTarget,:) ...
        +GkPrime(sourceIndex,:)+GlPrime(sourceIndex,:);
    dGl(range,:) = GlDouble(retainedTarget,:)-GlPrime(sourceIndex,:);
    ptr = ptr+number;
end
source = source(1:ptr);
target = target(1:ptr);
baseQ = baseQ(1:ptr,:);
dGl = dGl(1:ptr,:);

dCoeff = reciprocalCoefficients(dGl,stack.B{l});
if isfinite(options.spectatorModeCutoff)
    retained = max(abs(dCoeff),[],2) <= options.spectatorModeCutoff;
    source = source(retained);
    target = target(retained);
    baseQ = baseQ(retained,:);
    dGl = dGl(retained,:);
    dCoeff = dCoeff(retained,:);
end

data.layers = [j,k];
data.spectatorLayer = l;
data.dofIndexJ = indexJ;
data.dofIndexK = indexK;
data.numStates = [nJ,nK];
data.source = source;
data.target = target;
data.baseQ = baseQ;
data.dGl = dGl;
data.GkPrime = GkPrime;
data.GjDouble = GjDouble;
data.stackTau = stack.tau;
data.normalization = ...
    sqrt(abs(det(stack.B{j}))*abs(det(stack.B{k})))/(2*pi)^2;

if isempty(source)
    data.modeIndex = zeros(0,1);
    data.uniqueDGl = zeros(0,2);
    data.weightedHoppingByMode = cell(2,2);
    data.x = zeros(2,0);
    data.configurationGridSize = 0;
    data.recommendedConfigurationGridSize = 0;
    return
end

[uniqueCoeff,~,modeIndex] = unique(dCoeff,'rows','stable');
uniqueDGl = (stack.B{l}*uniqueCoeff.').';
maximumMode = max(abs(uniqueCoeff),[],'all');
recommendedGridSize = max( ...
    options.minimumConfigurationGridSize,2*maximumMode+1);
if isempty(options.configurationGridSize)
    configurationGridSize = recommendedGridSize;
else
    configurationGridSize = options.configurationGridSize;
    if configurationGridSize < recommendedGridSize
        warning('prepareRelaxedInterlayerHamiltonian:configurationAliasing', ...
            ['configurationGridSize=%d is below the recommended size %d ' ...
             'for pair (%d,%d). Check convergence carefully.'], ...
            configurationGridSize,recommendedGridSize,j,k);
    end
end

[x,xWeight] = polarQuadrature( ...
    options.realSpaceRadialOrder,options.realSpaceAngularOrder, ...
    options.realSpaceCutoff);
[bCell,bWeight] = periodicCellTrapezoid( ...
    stack.A{l},configurationGridSize);

nX = size(x,2);
nModes = size(uniqueDGl,1);
weightedHoppingByMode = cell(2,2);
for alpha = 1:2
    for beta = 1:2
        relaxedPosition = evaluateRelaxedInterlayerPosition( ...
            stack,fields,x,bCell,j,k,l,alpha,beta,options.useGPU);
        h = realSpaceInterlayerHopping( ...
            relaxedPosition,stack,j,k,alpha,beta);
        hoppingSamples = reshape(h,nX,size(bCell,2));
        coefficient = complex(zeros(nX,nModes));
        hoppingDevice = toDevice(hoppingSamples,options.useGPU);
        bDevice = toDevice(bCell,options.useGPU);
        bWeightDevice = toDevice(bWeight(:).',options.useGPU);
        for firstMode = 1:options.modeBatchSize:nModes
            lastMode = min(firstMode+options.modeBatchSize-1,nModes);
            modeRange = firstMode:lastMode;
            dGlDevice = toDevice(uniqueDGl(modeRange,:),options.useGPU);
            phaseB = exp(1i*(dGlDevice*bDevice)).*bWeightDevice;
            coefficient(:,modeRange) = toHost( ...
                hoppingDevice*phaseB.',options.useGPU);
        end
        weightedHoppingByMode{alpha,beta} = xWeight(:).*coefficient;
        clear relaxedPosition h hoppingSamples hoppingDevice coefficient
    end
end

data.modeIndex = modeIndex;
data.uniqueDGl = uniqueDGl;
data.weightedHoppingByMode = weightedHoppingByMode;
data.x = x;
data.configurationGridSize = configurationGridSize;
data.recommendedConfigurationGridSize = recommendedGridSize;

if options.verbose
    fprintf(['Prepared interlayer (%d,%d) cache: %d candidate pairs, ' ...
        '%d spectator modes, %d x %d cell grid.\n'], ...
        j,k,numel(source),nModes,configurationGridSize,configurationGridSize);
end
end

function options = mergePreservingEmpty(options,defaults,emptyName)
if nargin < 1 || isempty(options)
    options = struct;
end
names = fieldnames(defaults);
for n = 1:numel(names)
    name = names{n};
    if ~isfield(options,name) || ...
            (isempty(options.(name)) && ~strcmp(name,emptyName))
        options.(name) = defaults.(name);
    end
end
end

function validateOptions(options)
integerNames = {'realSpaceRadialOrder','realSpaceAngularOrder', ...
    'minimumConfigurationGridSize','modeBatchSize','pairBatchSize'};
for n = 1:numel(integerNames)
    value = options.(integerNames{n});
    assert(value >= 1 && value == fix(value), ...
        'options.%s must be a positive integer.',integerNames{n});
end
if ~isempty(options.configurationGridSize)
    assert(options.configurationGridSize >= 1 && ...
        options.configurationGridSize == fix(options.configurationGridSize), ...
        'configurationGridSize must be empty or a positive integer.');
end
assert(options.realSpaceCutoff > 0,'realSpaceCutoff must be positive.');
assert(options.spectatorModeCutoff >= 0, ...
    'spectatorModeCutoff must be nonnegative.');
assert(options.momentumInnerRadius >= 0 && ...
    options.momentumOuterRadius > options.momentumInnerRadius, ...
    'Momentum radii must satisfy 0 <= inner < outer.');
end
