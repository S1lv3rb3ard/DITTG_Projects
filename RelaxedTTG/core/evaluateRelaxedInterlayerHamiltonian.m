function [Hinter,blocks,info] = evaluateRelaxedInterlayerHamiltonian(cache,q)
%EVALUATERELAXEDINTERLAYERHAMILTONIAN Evaluate cached interlayer terms.
assert(isequal(size(q),[2,1]),'q must be a 2-by-1 column vector.');
distance = norm(q-cache.qReference);
assert(distance <= cache.qRadius+1e-12*max(1,cache.qRadius), ...
    ['q lies outside the domain prepared in the cache. Rebuild the cache ' ...
     'with all intended DoS/LDoS q-points or a larger qPadding.']);

upper = sparse(cache.dimension,cache.dimension);
blocks = cell(3,3);
pairInfo = cell(1,2);
for pair = 1:2
    data = cache.pair{pair};
    [block,pairInfo{pair}] = evaluatePair(data,q,cache.options);
    j = data.layers(1);
    k = data.layers(2);
    blocks{j,k} = block;
    blocks{k,j} = block';
    orbitalJ = reshape( ...
        [2*data.dofIndexJ-1,2*data.dofIndexJ].',[],1);
    orbitalK = reshape( ...
        [2*data.dofIndexK-1,2*data.dofIndexK].',[],1);
    upper(orbitalJ,orbitalK) = block;
end

Hinter = upper+upper';
info.pair = pairInfo;
info.numNonzeros = nnz(Hinter);
info.relativeHermiticityError = ...
    norm(Hinter-Hinter','fro')/max(norm(Hinter,'fro'),eps);
end

function [block,info] = evaluatePair(data,q,options)
nJ = data.numStates(1);
nK = data.numStates(2);
Qall = data.baseQ+q.';
chi = smoothCutoff(vecnorm(Qall,2,2), ...
    options.momentumInnerRadius,options.momentumOuterRadius);
retained = find(chi > 0);
if isempty(retained)
    block = sparse(2*nJ,2*nK);
    info = makeInfo(data,0,block);
    return
end

source = data.source(retained);
target = data.target(retained);
Qall = Qall(retained,:);
dGl = data.dGl(retained,:);
modeIndex = data.modeIndex(retained);
chi = chi(retained);
number = numel(retained);
valuesByOrbital = complex(zeros(number,4));

j = data.layers(1);
k = data.layers(2);
column = 0;
for alpha = 1:2
    for beta = 1:2
        column = column+1;
        integral = complex(zeros(number,1));
        modeList = unique(modeIndex,'stable').';
        weightedHopping = data.weightedHoppingByMode{alpha,beta};
        xDevice = toDevice(data.x,options.useGPU);
        for mode = modeList
            localPairs = find(modeIndex == mode);
            hDevice = toDevice(weightedHopping(:,mode),options.useGPU);
            for firstPair = 1:options.pairBatchSize:numel(localPairs)
                lastPair = min(firstPair+options.pairBatchSize-1, ...
                    numel(localPairs));
                pairRange = localPairs(firstPair:lastPair);
                QDevice = toDevice(Qall(pairRange,:),options.useGPU);
                integral(pairRange) = toHost( ...
                    exp(-1i*(QDevice*xDevice))*hDevice,options.useGPU);
            end
        end

        tauJ = (alpha-1)*dataTau(data,j);
        tauK = (beta-1)*dataTau(data,k);
        phase = exp(1i*(dataGjDouble(data,target)*tauJ)) .* ...
            exp(1i*((dGl-dataGkPrime(data,source))*tauK));
        valuesByOrbital(:,column) = ...
            data.normalization*chi.*phase.*integral;
    end
end

rowIndex = zeros(4*number,1);
colIndex = zeros(4*number,1);
flatValues = complex(zeros(4*number,1));
rowBase = 2*(source-1);
colBase = 2*(target-1);
rowIndex(1:4:end) = rowBase+1;
rowIndex(2:4:end) = rowBase+1;
rowIndex(3:4:end) = rowBase+2;
rowIndex(4:4:end) = rowBase+2;
colIndex(1:4:end) = colBase+1;
colIndex(2:4:end) = colBase+2;
colIndex(3:4:end) = colBase+1;
colIndex(4:4:end) = colBase+2;
flatValues(1:4:end) = valuesByOrbital(:,1);
flatValues(2:4:end) = valuesByOrbital(:,2);
flatValues(3:4:end) = valuesByOrbital(:,3);
flatValues(4:4:end) = valuesByOrbital(:,4);
block = sparse(rowIndex,colIndex,flatValues,2*nJ,2*nK);
info = makeInfo(data,number,block);
end

function value = dataTau(data,layer)
value = data.stackTau(:,layer);
end

function value = dataGjDouble(data,index)
value = data.GjDouble(index,:);
end

function value = dataGkPrime(data,index)
value = data.GkPrime(index,:);
end

function info = makeInfo(data,number,block)
info.layers = data.layers;
info.spectatorLayer = data.spectatorLayer;
info.dofIndexJ = data.dofIndexJ;
info.dofIndexK = data.dofIndexK;
info.numStates = data.numStates;
info.numRetainedStatePairs = number;
info.numCandidateStatePairs = numel(data.source);
info.numUniqueSpectatorModes = size(data.uniqueDGl,1);
info.configurationGridSize = data.configurationGridSize;
info.recommendedConfigurationGridSize = ...
    data.recommendedConfigurationGridSize;
info.numNonzeros = nnz(block);
end
