function relaxedPosition = evaluateRelaxedInterlayerPosition( ...
    stack,fields,x,bCell,j,k,l,alpha,beta,useGPU)
%EVALUATERELAXEDINTERLAYERPOSITION Construct x+u_j-u_k on x-by-b nodes.
%
% For (j,k,l)=(1,2,3), this evaluates
% x + u_1(x-tau_2beta,x+tau_2beta+b_3)
%   - u_2(-x+tau_1alpha,b_3).
if nargin < 10
    useGPU = false;
end

nX = size(x,2);
nB = size(bCell,2);
xDevice = toDevice(x,useGPU);
bDevice = toDevice(bCell,useGPU);
xAll = repmat(xDevice,1,nB);
bAll = repelem(bDevice,1,nX);
tauJ = toDevice((alpha-1)*stack.tau(:,j),useGPU);
tauK = toDevice((beta-1)*stack.tau(:,k),useGPU);

coordinatesJ = cell(1,3);
coordinatesJ{k} = xAll-tauK;
coordinatesJ{l} = xAll+tauK+bAll;
uJ = evaluateField(fields{j},j,coordinatesJ);

coordinatesK = cell(1,3);
coordinatesK{j} = -xAll+tauJ;
coordinatesK{l} = bAll;
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
