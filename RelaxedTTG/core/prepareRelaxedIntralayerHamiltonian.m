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

shellMap = stack.A{layer}/stack.A{options.shellReferenceLayer};
RA = shellMap*shellsReference.matA;
RB = shellMap*shellsReference.matB;
t = stack.tau(:,layer);

N = options.configurationGridSize;
[bk1,bk2,bl1,bl2] = ndgrid((0:N-1)/N);
bk = stack.A{kLayer}*[bk1(:),bk2(:)].';
bl = stack.A{lLayer}*[bl1(:),bl2(:)].';
clear bk1 bk2 bl1 bl2

useGPU = options.useGPU;
bkDevice = toDevice(bk,useGPU);
blDevice = toDevice(bl,useGPU);
RADevice = toDevice(RA,useGPU);
RBDevice = toDevice(RB,useGPU);
tDevice = toDevice(t,useGPU);
weight = 1/N^4;

RAcell = num2cell(RADevice,1);
RBcell = num2cell(RBDevice,1);
inputAA = cellfun(@(R) ...
    R+relaxation(R+bkDevice,R+blDevice)-relaxation(bkDevice,blDevice), ...
    RAcell,'UniformOutput',false);
inputAA = squeeze(vecnorm(cat(3,inputAA{:}),2,1));
hopAA = toDevice(shellsReference.intraAA(toHost(inputAA,useGPU)),useGPU);

inputBB = cellfun(@(R) ...
    R+relaxation(R+bkDevice+tDevice,R+blDevice+tDevice) ...
     -relaxation(bkDevice+tDevice,blDevice+tDevice), ...
    RAcell,'UniformOutput',false);
inputBB = squeeze(vecnorm(cat(3,inputBB{:}),2,1));
hopBB = toDevice(shellsReference.intraAA(toHost(inputBB,useGPU)),useGPU);

inputAB = cellfun(@(R) ...
    R+relaxation(R+bkDevice,R+blDevice) ...
     -relaxation(bkDevice+tDevice,blDevice+tDevice), ...
    RBcell,'UniformOutput',false);
inputAB = squeeze(vecnorm(cat(3,inputAB{:}),2,1));
hopAB = toDevice(shellsReference.intraAB(toHost(inputAB,useGPU)),useGPU);

inputBA = cellfun(@(R) ...
    R+relaxation(R+bkDevice+tDevice,R+blDevice+tDevice) ...
     -relaxation(bkDevice,blDevice), ...
    RBcell,'UniformOutput',false);
inputBA = squeeze(vecnorm(cat(3,inputBA{:}),2,1));
hopBA = toDevice(shellsReference.intraAB(toHost(inputBA,useGPU)),useGPU);
clear inputAA inputBB inputAB inputBA RAcell RBcell

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
coeffAA = complex(zeros(pairCapacity,size(RA,2)));
coeffBB = complex(zeros(pairCapacity,size(RA,2)));
coeffAB = complex(zeros(pairCapacity,size(RB,2)));
coeffBA = complex(zeros(pairCapacity,size(RB,2)));
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
    baseQ = Gk(sourceIndex,:)+Gl(sourceIndex,:);
    baseQ2 = baseQ+dGk+dGl;

    dGkDevice = toDevice(dGk,useGPU);
    dGlDevice = toDevice(dGl,useGPU);
    baseQDevice = toDevice(baseQ,useGPU);
    baseQ2Device = toDevice(baseQ2,useGPU);
    phaseG = exp(1i*dGkDevice*bkDevice) ...
        .*exp(1i*dGlDevice*blDevice);
    phaseRA = exp(-1i*baseQDevice*RADevice);
    phaseRB = exp(-1i*baseQDevice*RBDevice);

    coeffAA(range,:) = toHost( ...
        weight*phaseRA.*(phaseG*hopAA),useGPU);
    coeffBB(range,:) = toHost( ...
        weight*exp(1i*(dGkDevice+dGlDevice)*tDevice) ...
        .*phaseRA.*(phaseG*hopBB),useGPU);
    coeffAB(range,:) = toHost( ...
        weight*exp(1i*baseQ2Device*tDevice) ...
        .*phaseRB.*(phaseG*hopAB),useGPU);
    coeffBA(range,:) = toHost( ...
        weight*exp(-1i*baseQDevice*tDevice) ...
        .*phaseRB.*(phaseG*hopBA),useGPU);
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
data.RA = RA;
data.RBminusTau = RB-t;
data.RBplusTau = RB+t;
data.numRetainedUnorderedPairs = ptr;

if options.verbose
    fprintf('Prepared intralayer %d cache: %d unordered pairs.\n', ...
        layer,ptr);
end
end
