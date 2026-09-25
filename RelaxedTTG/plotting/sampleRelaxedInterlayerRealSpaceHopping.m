function data = sampleRelaxedInterlayerRealSpaceHopping( ...
    stack,relaxationFields,pair,options)
%SAMPLERELAXEDINTERLAYERREALSPACEHOPPING Sample all four orbital channels.
if nargin < 4
    options = struct;
end
defaults = struct( ...
    'extent',4*stack.aG, ...
    'gridSize',121, ...
    'spectatorFraction',[0;0], ...
    'useGPU',false);
options = mergeOptions(options,defaults);

[j,k,l] = validatePair(pair);
assert(options.gridSize >= 2 && options.gridSize == fix(options.gridSize), ...
    'gridSize must be an integer at least two.');
fraction = options.spectatorFraction(:);
assert(numel(fraction) == 2, ...
    'spectatorFraction must contain two components.');

axisValues = linspace(-options.extent,options.extent,options.gridSize);
[x1,x2] = ndgrid(axisValues,axisValues);
x = [x1(:),x2(:)].';
b = stack.A{l}*fraction;

data.layers = [j,k];
data.spectatorLayer = l;
data.spectatorFraction = fraction;
data.axis = axisValues;
data.x1 = x1;
data.x2 = x2;
data.hopping = cell(2,2);
data.relaxedPosition = cell(2,2);

for alpha = 1:2
    for beta = 1:2
        relaxedPosition = evaluateRelaxedInterlayerPosition( ...
            stack,relaxationFields,x,b,j,k,l,alpha,beta,options.useGPU);
        hopping = realSpaceInterlayerHopping( ...
            relaxedPosition,stack,j,k,alpha,beta);
        data.relaxedPosition{alpha,beta} = reshape( ...
            toHost(relaxedPosition,options.useGPU),2,options.gridSize, ...
            options.gridSize);
        data.hopping{alpha,beta} = reshape( ...
            toHost(hopping,options.useGPU),options.gridSize,options.gridSize);
    end
end
end

function [j,k,l] = validatePair(pair)
pair = pair(:).';
assert(isequal(pair,[1,2]) || isequal(pair,[2,3]), ...
    'pair must be [1,2] or [2,3].');
j = pair(1);
k = pair(2);
l = setdiff(1:3,pair);
end
