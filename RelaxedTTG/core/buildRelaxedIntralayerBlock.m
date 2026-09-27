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

N = options.configurationGridSize;
[bk1,bk2,bl1,bl2] = ndgrid((0:N-1)/N);
bk = stack.A{kLayer}*[bk1(:),bk2(:)].';
bl = stack.A{lLayer}*[bl1(:),bl2(:)].';
clear bk1 bk2 bl1 bl2

useGPU = options.useGPU;
bk = toDevice(bk,useGPU);
bl = toDevice(bl,useGPU);
t = stack.tau(:,layer);
tDevice = toDevice(t,useGPU);
weight = 1/N^4;

% Relaxed hoppings for AA, AB, BA, BB (row, column) at disregistry-consistent
% orbital configurations; bond0 = d + tau_row - tau_col.
channels = sampleRelaxedIntralayerChannels(stack,shellsReference, ...
    relaxation,layer,bk,bl, ...
    struct('shellReferenceLayer',options.shellReferenceLayer, ...
           'useGPU',useGPU));
bondAADevice = toDevice(channels(1).bond0,useGPU);
bondABDevice = toDevice(channels(2).bond0,useGPU);
bondBADevice = toDevice(channels(3).bond0,useGPU);
bondBBDevice = toDevice(channels(4).bond0,useGPU);
hopAA = channels(1).hopping;
hopAB = channels(2).hopping;
hopBA = channels(3).hopping;
hopBB = channels(4).hopping;
clear channels

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

    % Theorem 4.2 with Q = q + S G' and bond0 = d + tau_row - tau_col.
    Q = q.'+Gk(i,:)+Gl(i,:);
    dGkDevice = toDevice(dGk,useGPU);
    dGlDevice = toDevice(dGl,useGPU);
    QDevice = toDevice(Q,useGPU);

    phaseG = exp(1i*dGkDevice*bk).*exp(1i*dGlDevice*bl);
    phaseTauB = exp(1i*(dGkDevice+dGlDevice)*tDevice);

    hAA = weight*sum(exp(-1i*QDevice*bondAADevice) ...
        .*(phaseG*hopAA),2);
    hAB = weight*phaseTauB.*sum(exp(-1i*QDevice*bondABDevice) ...
        .*(phaseG*hopAB),2);
    hBA = weight*sum(exp(-1i*QDevice*bondBADevice) ...
        .*(phaseG*hopBA),2);
    hBB = weight*phaseTauB.*sum(exp(-1i*QDevice*bondBBDevice) ...
        .*(phaseG*hopBB),2);
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
