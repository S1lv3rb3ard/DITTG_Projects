function [Hinter,blocks,info] = buildRelaxedInterlayerHamiltonian( ...
    stack,DoF,q,relaxationFields,options)
%BUILDRELAXEDINTERLAYERHAMILTONIAN Assemble relaxed (1,2) and (2,3) blocks.
%
% The R^2 Fourier integral is evaluated by direct polar quadrature.  The
% periodic spectator-cell integral is evaluated by a direct trapezoidal
% sum.  No FFT is used.  State pairs are grouped by G_l''-G_l' so each
% required spectator-cell coefficient is integrated once and reused.
if nargin < 5
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
options = mergeOptionsPreserveEmpty(options,defaults,'configurationGridSize');
validateOptions(options);

assert(isequal(size(q),[2,1]),'q must be a 2-by-1 column vector.');
assert(iscell(relaxationFields) && numel(relaxationFields) == 3, ...
    'relaxationFields must contain three function handles.');
for layer = 1:3
    assert(isa(relaxationFields{layer},'function_handle'), ...
        'relaxationFields{%d} must be a function handle.',layer);
end

nDoF = size(DoF,1);
upper = sparse(2*nDoF,2*nDoF);
blocks = cell(3,3);
pairList = [1,2,3;2,3,1];
pairInfo = cell(1,2);

for pair = 1:2
    j = pairList(pair,1);
    k = pairList(pair,2);
    l = pairList(pair,3);
    [block,pairInfo{pair}] = buildPairBlock( ...
        stack,DoF,q,relaxationFields,options,j,k,l);
    blocks{j,k} = block;
    blocks{k,j} = block';

    indexJ = pairInfo{pair}.dofIndexJ;
    indexK = pairInfo{pair}.dofIndexK;
    orbitalJ = reshape([2*indexJ-1,2*indexJ].',[],1);
    orbitalK = reshape([2*indexK-1,2*indexK].',[],1);
    upper(orbitalJ,orbitalK) = block;
end

Hinter = upper+upper';
relativeError = norm(Hinter-Hinter','fro')/max(norm(Hinter,'fro'),eps);
assert(relativeError < 1e-13, ...
    'The assembled interlayer Hamiltonian is not Hermitian.');
info.pair = pairInfo;
info.options = options;
info.numNonzeros = nnz(Hinter);
info.relativeHermiticityError = relativeError;
end

function [block,info] = buildPairBlock( ...
    stack,DoF,q,fields,options,j,k,l)
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

% Enumerate only pairs inside the Q cutoff.
pairCount = 0;
for sourceIndex = 1:nJ
    Q = q.'+GjDouble+GkPrime(sourceIndex,:)+GlPrime(sourceIndex,:);
    pairCount = pairCount+nnz(vecnorm(Q,2,2) < options.momentumOuterRadius);
end

source = zeros(pairCount,1);
target = zeros(pairCount,1);
Qall = zeros(pairCount,2);
dGl = zeros(pairCount,2);
chi = zeros(pairCount,1);
ptr = 0;
for sourceIndex = 1:nJ
    Q = q.'+GjDouble+GkPrime(sourceIndex,:)+GlPrime(sourceIndex,:);
    normQ = vecnorm(Q,2,2);
    retainedTarget = find(normQ < options.momentumOuterRadius);
    number = numel(retainedTarget);
    if number == 0
        continue
    end
    range = ptr+(1:number);
    source(range) = sourceIndex;
    target(range) = retainedTarget;
    Qall(range,:) = Q(retainedTarget,:);
    dGl(range,:) = GlDouble(retainedTarget,:)-GlPrime(sourceIndex,:);
    chi(range) = smoothCutoff(normQ(retainedTarget), ...
        options.momentumInnerRadius,options.momentumOuterRadius);
    ptr = ptr+number;
end

retained = chi > 0;
source = source(retained);
target = target(retained);
Qall = Qall(retained,:);
dGl = dGl(retained,:);
chi = chi(retained);

dCoeff = reciprocalCoefficients(dGl,stack.B{l});
if isfinite(options.spectatorModeCutoff)
    retained = max(abs(dCoeff),[],2) <= options.spectatorModeCutoff;
    source = source(retained);
    target = target(retained);
    Qall = Qall(retained,:);
    dGl = dGl(retained,:);
    dCoeff = dCoeff(retained,:);
    chi = chi(retained);
end
pairCount = numel(chi);

if pairCount == 0
    block = sparse(2*nJ,2*nK);
    info = makePairInfo(j,k,l,indexJ,indexK,nJ,nK,0,0,0,0,block);
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
        warning('buildRelaxedInterlayerHamiltonian:configurationAliasing', ...
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

blockValues = complex(zeros(pairCount,4));
normalization = sqrt(abs(det(stack.B{j}))*abs(det(stack.B{k})))/(2*pi)^2;
modePairs = accumarray(modeIndex,(1:pairCount)', ...
    [size(uniqueCoeff,1),1],@(indices){indices});

column = 0;
for alpha = 1:2
    for beta = 1:2
        column = column+1;
        hoppingSamples = configurationHoppingSamples( ...
            stack,fields,x,bCell,j,k,l,alpha,beta,options.useGPU);
        integral = evaluateGroupedFourierIntegral( ...
            hoppingSamples,x,xWeight,bCell,bWeight,Qall, ...
            uniqueDGl,modePairs,options);

        tauJ = (alpha-1)*stack.tau(:,j);
        tauK = (beta-1)*stack.tau(:,k);
        phase = exp(1i*(GjDouble(target,:)*tauJ)) .* ...
            exp(1i*((dGl-GkPrime(source,:))*tauK));
        blockValues(:,column) = normalization*chi.*phase.*integral;
        clear hoppingSamples integral
    end
end

rowIndex = zeros(4*pairCount,1);
colIndex = zeros(4*pairCount,1);
values = complex(zeros(4*pairCount,1));
for entry = 1:pairCount
    rows = 2*(source(entry)-1)+(1:2);
    columns = 2*(target(entry)-1)+(1:2);
    range = 4*(entry-1)+(1:4);
    rowIndex(range) = [rows(1);rows(1);rows(2);rows(2)];
    colIndex(range) = [columns(1);columns(2);columns(1);columns(2)];
    values(range) = blockValues(entry,:).';
end
block = sparse(rowIndex,colIndex,values,2*nJ,2*nK);

info = makePairInfo(j,k,l,indexJ,indexK,nJ,nK,pairCount, ...
    size(uniqueCoeff,1),configurationGridSize,recommendedGridSize,block);
if options.verbose
    fprintf(['Interlayer (%d,%d): %d retained pairs, %d unique spectator ' ...
        'modes, %d x %d cell grid, nnz=%d.\n'], ...
        j,k,pairCount,size(uniqueCoeff,1),configurationGridSize, ...
        configurationGridSize,nnz(block));
end
end

function hoppingSamples = configurationHoppingSamples( ...
    stack,fields,x,bCell,j,k,l,alpha,beta,useGPU)
nX = size(x,2);
nB = size(bCell,2);
relaxedPosition = evaluateRelaxedInterlayerPosition( ...
    stack,fields,x,bCell,j,k,l,alpha,beta,useGPU);
h = realSpaceInterlayerHopping( ...
    relaxedPosition,stack,j,k,alpha,beta);
hoppingSamples = reshape(h,nX,nB);
end

function integral = evaluateGroupedFourierIntegral( ...
    hoppingSamples,x,xWeight,bCell,bWeight,Qall,uniqueDGl,modePairs,options)
% First integrate each required spectator mode over Gamma_l, then reuse it
% for every state pair with that mode.  Both stages are direct quadrature.
nPairs = size(Qall,1);
nModes = size(uniqueDGl,1);
integral = complex(zeros(nPairs,1));
useGPU = options.useGPU;
xDevice = toDevice(x,useGPU);
xWeightDevice = toDevice(xWeight(:).',useGPU);
bDevice = toDevice(bCell,useGPU);
bWeightDevice = toDevice(bWeight(:).',useGPU);

for firstMode = 1:options.modeBatchSize:nModes
    lastMode = min(firstMode+options.modeBatchSize-1,nModes);
    modeRange = firstMode:lastMode;
    dGlDevice = toDevice(uniqueDGl(modeRange,:),useGPU);
    phaseB = exp(1i*(dGlDevice*bDevice)).*bWeightDevice;
    hoppingByMode = hoppingSamples*phaseB.';

    for localMode = 1:numel(modeRange)
        globalMode = modeRange(localMode);
        pairList = modePairs{globalMode};
        for firstPair = 1:options.pairBatchSize:numel(pairList)
            lastPair = min(firstPair+options.pairBatchSize-1,numel(pairList));
            pairRange = pairList(firstPair:lastPair);
            QDevice = toDevice(Qall(pairRange,:),useGPU);
            phaseX = exp(-1i*(QDevice*xDevice)).*xWeightDevice;
            value = phaseX*hoppingByMode(:,localMode);
            integral(pairRange) = toHost(value,useGPU);
        end
    end
    clear phaseB hoppingByMode
end
end

function info = makePairInfo(j,k,l,indexJ,indexK,nJ,nK,pairCount, ...
    modeCount,gridSize,recommendedGridSize,block)
info.layers = [j,k];
info.spectatorLayer = l;
info.dofIndexJ = indexJ;
info.dofIndexK = indexK;
info.numStates = [nJ,nK];
info.numRetainedStatePairs = pairCount;
info.numUniqueSpectatorModes = modeCount;
info.configurationGridSize = gridSize;
info.recommendedConfigurationGridSize = recommendedGridSize;
info.numNonzeros = nnz(block);
end

function options = mergeOptionsPreserveEmpty(options,defaults,emptyName)
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
