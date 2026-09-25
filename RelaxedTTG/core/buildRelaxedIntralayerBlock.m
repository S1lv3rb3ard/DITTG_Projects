function [Hj,info] = buildRelaxedIntralayerBlock( ...
    stack,DoF,layer,q,shellsReference,relaxation,options)
%BUILDRELAXEDINTRALAYERBLOCK Relaxed reciprocal intralayer block for layer j.
defaults = struct( ...
    'configurationGridSize',5, ...
    'couplingInnerRadius',6.25, ...
    'couplingOuterRadius',6.50, ...
    'shellReferenceLayer',2, ...
    'useGPU',false, ...
    'verbose',true);
options = mergeOptions(options,defaults);

assert(ismember(layer,1:3) && layer == fix(layer), ...
    'layer must be 1, 2, or 3.');
assert(isequal(size(q),[2,1]),'q must be a 2-by-1 column vector.');
assert(isa(relaxation,'function_handle'), ...
    'relaxation must be a function handle.');
assert(options.configurationGridSize >= 1 && ...
    options.configurationGridSize == fix(options.configurationGridSize), ...
    'configurationGridSize must be a positive integer.');
assert(ismember(options.shellReferenceLayer,1:3), ...
    'shellReferenceLayer must be 1, 2, or 3.');

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
bk = toDevice(bk,useGPU);
bl = toDevice(bl,useGPU);
RADevice = toDevice(RA,useGPU);
RBDevice = toDevice(RB,useGPU);
tDevice = toDevice(t,useGPU);
weight = 1/N^4;

RAcell = num2cell(RADevice,1);
RBcell = num2cell(RBDevice,1);

inputAA = cellfun(@(R) ...
    R+relaxation(R+bk,R+bl)-relaxation(bk,bl), ...
    RAcell,'UniformOutput',false);
inputAA = squeeze(vecnorm(cat(3,inputAA{:}),2,1));
hopAA = toDevice(shellsReference.intraAA(toHost(inputAA,useGPU)),useGPU);

inputBB = cellfun(@(R) ...
    R+relaxation(R+bk+tDevice,R+bl+tDevice) ...
     -relaxation(bk+tDevice,bl+tDevice), ...
    RAcell,'UniformOutput',false);
inputBB = squeeze(vecnorm(cat(3,inputBB{:}),2,1));
hopBB = toDevice(shellsReference.intraAA(toHost(inputBB,useGPU)),useGPU);

inputAB = cellfun(@(R) ...
    R+relaxation(R+bk,R+bl) ...
     -relaxation(bk+tDevice,bl+tDevice), ...
    RBcell,'UniformOutput',false);
inputAB = squeeze(vecnorm(cat(3,inputAB{:}),2,1));
hopAB = toDevice(shellsReference.intraAB(toHost(inputAB,useGPU)),useGPU);

inputBA = cellfun(@(R) ...
    R+relaxation(R+bk+tDevice,R+bl+tDevice) ...
     -relaxation(bk,bl), ...
    RBcell,'UniformOutput',false);
inputBA = squeeze(vecnorm(cat(3,inputBA{:}),2,1));
hopBA = toDevice(shellsReference.intraAB(toHost(inputBA,useGPU)),useGPU);
clear inputAA inputBB inputAB inputBA RAcell RBcell

pairCount = 0;
for i = 1:m
    dGk = Gk(i:m,:)-Gk(i,:);
    dGl = Gl(i:m,:)-Gl(i,:);
    rho = sqrt(sum(dGk.^2,2)+sum(dGl.^2,2));
    pairCount = pairCount+nnz(rho < options.couplingOuterRadius);
end

rowIndex = zeros(4*pairCount,1);
colIndex = zeros(4*pairCount,1);
values = complex(zeros(4*pairCount,1));
ptr = 0;

for i = 1:m
    allJ = (i:m).';
    dGk = Gk(allJ,:)-Gk(i,:);
    dGl = Gl(allJ,:)-Gl(i,:);
    rho = sqrt(sum(dGk.^2,2)+sum(dGl.^2,2));
    keep = rho < options.couplingOuterRadius;
    J = allJ(keep);
    dGk = dGk(keep,:);
    dGl = dGl(keep,:);
    couplingChi = smoothCutoff(rho(keep), ...
        options.couplingInnerRadius,options.couplingOuterRadius);
    keep = couplingChi > 0;
    J = J(keep);
    dGk = dGk(keep,:);
    dGl = dGl(keep,:);
    couplingChi = couplingChi(keep);
    if isempty(J)
        continue
    end

    Q = q.'+Gk(i,:)+Gl(i,:);
    Q2 = Q+dGk+dGl;
    dGkDevice = toDevice(dGk,useGPU);
    dGlDevice = toDevice(dGl,useGPU);
    QDevice = toDevice(Q,useGPU);
    Q2Device = toDevice(Q2,useGPU);

    phaseG = exp(1i*dGkDevice*bk).*exp(1i*dGlDevice*bl);
    phaseQA = exp(-1i*QDevice*RADevice);
    phaseQB = exp(-1i*QDevice*RBDevice);

    hAA = weight*sum(phaseQA.*(phaseG*hopAA),2);
    hBB = weight*exp(1i*(dGkDevice+dGlDevice)*tDevice) ...
        .*sum(phaseQA.*(phaseG*hopBB),2);
    hAB = weight*exp(1i*Q2Device*tDevice) ...
        .*sum(phaseQB.*(phaseG*hopAB),2);
    hBA = weight*exp(-1i*QDevice*tDevice) ...
        .*sum(phaseQB.*(phaseG*hopBA),2);
    blockValues = couplingChi.*toHost([hAA,hAB,hBA,hBB],useGPU);

    for p = 1:numel(J)
        j = J(p);
        Bij = [blockValues(p,1),blockValues(p,2); ...
               blockValues(p,3),blockValues(p,4)];
        if i == j
            Bij = 0.25*(Bij+Bij');
        end

        rows = 2*(i-1)+(1:2);
        columns = 2*(j-1)+(1:2);
        index = ptr+(1:4);
        rowIndex(index) = [rows(1);rows(1);rows(2);rows(2)];
        colIndex(index) = [columns(1);columns(2);columns(1);columns(2)];
        values(index) = [Bij(1,1);Bij(1,2);Bij(2,1);Bij(2,2)];
        ptr = ptr+4;
    end

    if options.verbose && (mod(i,100) == 0 || i == m)
        fprintf('Intralayer %d: completed row %d of %d.\n',layer,i,m);
    end
end

upper = sparse(rowIndex(1:ptr),colIndex(1:ptr),values(1:ptr),2*m,2*m);
Hj = upper+upper';
info.layer = layer;
info.otherLayers = otherLayers;
info.dofIndex = dofIndex;
info.numStates = m;
info.numRetainedUnorderedPairs = ptr/4;
info.numNonzeros = nnz(Hj);
info.relativeHermiticityError = ...
    norm(Hj-Hj','fro')/max(norm(Hj,'fro'),eps);
end
