function cache = prepareRelaxedIntralayerHamiltonian( ...
    stack,DoF,shellsReference,relaxationFields,options)
%PREPARERELAXEDINTRALAYERHAMILTONIAN Cache relaxed intralayer coefficients.
if nargin < 5
    options = struct;
end
defaults = struct( ...
    'configurationGridSize',5, ...
    'couplingInnerRadius',6.25, ...
    'couplingOuterRadius',6.50, ...
    'shellReferenceLayer',2, ...
    'useGPU',false, ...
    'verbose',true);
options = mergeOptions(options,defaults);

assert(iscell(relaxationFields) && numel(relaxationFields) == 3, ...
    'relaxationFields must contain three function handles.');
assert(options.configurationGridSize >= 1 && ...
    options.configurationGridSize == fix(options.configurationGridSize), ...
    'configurationGridSize must be a positive integer.');

cache.layer = cell(1,3);
cache.options = options;
cache.dimension = 2*size(DoF,1);
for layer = 1:3
    assert(isa(relaxationFields{layer},'function_handle'), ...
        'relaxationFields{%d} must be a function handle.',layer);
    cache.layer{layer} = prepareLayer( ...
        stack,DoF,shellsReference,relaxationFields{layer},layer,options);
end
end

function data = prepareLayer( ...
    stack,DoF,shellsReference,relaxation,layer,options)
otherLayers = setdiff(1:3,layer,'stable');
kLayer = otherLayers(1);
lLayer = otherLayers(2);
dofIndex = find(DoF(:,9) == layer);
Gk = DoF(dofIndex,2*kLayer+(1:2));
Gl = DoF(dofIndex,2*lLayer+(1:2));
m = numel(dofIndex);
assert(m > 0,'Layer %d has no retained DoFs.',layer);

N = options.configurationGridSize;
[bk1,bk2,bl1,bl2] = ndgrid((0:N-1)/N);
bk = stack.A{kLayer}*[bk1(:),bk2(:)].';
bl = stack.A{lLayer}*[bl1(:),bl2(:)].';
clear bk1 bk2 bl1 bl2

useGPU = options.useGPU;
bkDevice = toDevice(bk,useGPU);
blDevice = toDevice(bl,useGPU);
weight = 1/N^4;
t = stack.tau(:,layer);
tDevice = toDevice(t,useGPU);

% Relaxed hoppings for the AA, AB, BA and BB channels (row, column), with
% every orbital displacement evaluated at its disregistry configuration.
channels = sampleRelaxedIntralayerChannels(stack,shellsReference, ...
    relaxation,layer,bkDevice,blDevice, ...
    struct('shellReferenceLayer',options.shellReferenceLayer, ...
           'useGPU',useGPU));
bondAA = channels(1).bond0;
bondAB = channels(2).bond0;
bondBA = channels(3).bond0;
bondBB = channels(4).bond0;
hopAA = channels(1).hopping;
hopAB = channels(2).hopping;
hopBA = channels(3).hopping;
hopBB = channels(4).hopping;
bondAADevice = toDevice(bondAA,useGPU);
bondABDevice = toDevice(bondAB,useGPU);
bondBADevice = toDevice(bondBA,useGPU);
bondBBDevice = toDevice(bondBB,useGPU);
clear channels

pairCapacity = 0;
for sourceIndex = 1:m
    dGk = Gk(sourceIndex:m,:)-Gk(sourceIndex,:);
    dGl = Gl(sourceIndex:m,:)-Gl(sourceIndex,:);
    rho = sqrt(sum(dGk.^2,2)+sum(dGl.^2,2));
    pairCapacity = pairCapacity+nnz(rho < options.couplingOuterRadius);
end

source = zeros(pairCapacity,1);
target = zeros(pairCapacity,1);
couplingChi = zeros(pairCapacity,1);
coeffAA = complex(zeros(pairCapacity,size(bondAA,2)));
coeffBB = complex(zeros(pairCapacity,size(bondBB,2)));
coeffAB = complex(zeros(pairCapacity,size(bondAB,2)));
coeffBA = complex(zeros(pairCapacity,size(bondBA,2)));
ptr = 0;

for sourceIndex = 1:m
    allTarget = (sourceIndex:m).';
    dGk = Gk(allTarget,:)-Gk(sourceIndex,:);
    dGl = Gl(allTarget,:)-Gl(sourceIndex,:);
    rho = sqrt(sum(dGk.^2,2)+sum(dGl.^2,2));
    chi = smoothCutoff(rho, ...
        options.couplingInnerRadius,options.couplingOuterRadius);
    retained = chi > 0;
    allTarget = allTarget(retained);
    dGk = dGk(retained,:);
    dGl = dGl(retained,:);
    chi = chi(retained);
    number = numel(allTarget);
    if number == 0
        continue
    end

    range = ptr+(1:number);
    source(range) = sourceIndex;
    target(range) = allTarget;
    couplingChi(range) = chi;
    % Theorem 4.2: element (G',j alpha),(G'',j beta) equals
    %   sum_d [h]_{G''-G'}(d) exp(-i(q+SG').bond0) exp(i S(G''-G').tau_beta),
    % with bond0 = d + tau_alpha - tau_beta. The q-dependent factor
    % exp(-i q.bond0) is applied in evaluateRelaxedIntralayerHamiltonian.
    baseQ = Gk(sourceIndex,:)+Gl(sourceIndex,:);

    dGkDevice = toDevice(dGk,useGPU);
    dGlDevice = toDevice(dGl,useGPU);
    baseQDevice = toDevice(baseQ,useGPU);
    phaseG = exp(1i*dGkDevice*bkDevice) ...
        .*exp(1i*dGlDevice*blDevice);
    phaseTauB = exp(1i*(dGkDevice+dGlDevice)*tDevice);

    coeffAA(range,:) = toHost(weight ...
        *exp(-1i*baseQDevice*bondAADevice).*(phaseG*hopAA),useGPU);
    coeffAB(range,:) = toHost(weight*phaseTauB ...
        .*exp(-1i*baseQDevice*bondABDevice).*(phaseG*hopAB),useGPU);
    coeffBA(range,:) = toHost(weight ...
        *exp(-1i*baseQDevice*bondBADevice).*(phaseG*hopBA),useGPU);
    coeffBB(range,:) = toHost(weight*phaseTauB ...
        .*exp(-1i*baseQDevice*bondBBDevice).*(phaseG*hopBB),useGPU);
    ptr = ptr+number;
end

source = source(1:ptr);
target = target(1:ptr);
data.layer = layer;
data.dofIndex = dofIndex;
data.numStates = m;
data.source = source;
data.target = target;
data.diagonal = source == target;
data.couplingChi = couplingChi(1:ptr);
data.coeffAA = coeffAA(1:ptr,:);
data.coeffBB = coeffBB(1:ptr,:);
data.coeffAB = coeffAB(1:ptr,:);
data.coeffBA = coeffBA(1:ptr,:);
data.bondAA = bondAA;
data.bondAB = bondAB;
data.bondBA = bondBA;
data.bondBB = bondBB;
data.numRetainedUnorderedPairs = ptr;

if options.verbose
    fprintf('Prepared intralayer %d cache: %d unordered pairs.\n', ...
        layer,ptr);
end
end
