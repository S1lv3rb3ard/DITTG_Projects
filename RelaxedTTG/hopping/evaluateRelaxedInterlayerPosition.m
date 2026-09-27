function relaxedPosition = evaluateRelaxedInterlayerPosition( ...
    stack,fields,x,bCell,j,k,l,alpha,beta,useGPU)
%EVALUATERELAXEDINTERLAYERPOSITION Construct x+u_{j,alpha}-u_{k,beta} on x-by-b nodes.
%
% x is the continuous interlayer coordinate d_tau = R - R' + tau_{j alpha}
% - tau_{k beta} of the bond from (R',k beta) to (R,j alpha), and bCell is
% the spectator configuration of the column lattice point, b = R' mod R_l.
% With e = x - tau_{j alpha} + tau_{k beta} = R - R', the lattice points have
% configurations T_j R = (e, e + b) in slots (k,l) and T_k R' = (-e, b) in
% slots (j,l). Each orbital field is evaluated at its disregistry
% configuration u_{j,alpha}(c) = u_j(c + D_j tau_{j alpha}):
%
%   u_{j,alpha}: slot k  e + (I - A_k A_j^{-1}) tau_{j alpha}
%                slot l  e + b + (I - A_l A_j^{-1}) tau_{j alpha}
%   u_{k,beta} : slot j -e + (I - A_j A_k^{-1}) tau_{k beta}
%                slot l  b + (I - A_l A_k^{-1}) tau_{k beta}
%
% See getDisregistryOffsets and Section 3.2 of the accompanying paper.
if nargin < 10
    useGPU = false;
end

nX = size(x,2);
nB = size(bCell,2);
xDevice = toDevice(x,useGPU);
bDevice = toDevice(bCell,useGPU);
xAll = repmat(xDevice,1,nB);
bAll = repelem(bDevice,1,nX);
tauJ = (alpha-1)*stack.tau(:,j);
tauK = (beta-1)*stack.tau(:,k);
offsetJ = toDevice(getDisregistryOffsets(stack,j,tauJ),useGPU);
offsetK = toDevice(getDisregistryOffsets(stack,k,tauK),useGPU);
eAll = xAll-toDevice(tauJ,useGPU)+toDevice(tauK,useGPU);

coordinatesJ = cell(1,3);
coordinatesJ{k} = eAll+offsetJ(:,k);
coordinatesJ{l} = eAll+bAll+offsetJ(:,l);
uJ = evaluateField(fields{j},j,coordinatesJ);

coordinatesK = cell(1,3);
coordinatesK{j} = -eAll+offsetK(:,j);
coordinatesK{l} = bAll+offsetK(:,l);
uK = evaluateField(fields{k},k,coordinatesK);
relaxedPosition = xAll+uJ-uK;
end

function u = evaluateField(field,layer,coordinates)
otherLayers = setdiff(1:3,layer,'stable');
assert(~isempty(coordinates{otherLayers(1)}) && ...
       ~isempty(coordinates{otherLayers(2)}), ...
    'Missing configuration coordinate for layer %d.',layer);
u = field(coordinates{otherLayers(1)},coordinates{otherLayers(2)});
assert(isequal(size(u),size(coordinates{otherLayers(1)})), ...
    'relaxationFields{%d} must return a 2-by-N array.',layer);
end
