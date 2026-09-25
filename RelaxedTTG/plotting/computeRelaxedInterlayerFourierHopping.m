function values = computeRelaxedInterlayerFourierHopping( ...
    stack,relaxationFields,pair,alpha,beta,Q,spectatorMode,options)
%COMPUTERELAXEDINTERLAYERFOURIERHOPPING Direct continuous Fourier integral.
%
% transformType='spectatorMode' and dG_l=B_l*spectatorMode returns
% (2*pi)^(-2) avg_{Gamma_l} int_R2 h_rel(x,b_l)
% exp(-i*Q*x) exp(i*dG_l*b_l) dx db_l.
%
% transformType='spectatorSlice' instead fixes
% b_l=A_l*spectatorFraction and transforms that two-dimensional slice:
% (2*pi)^(-2) int_R2 h_rel(x,b_l) exp(-i*Q*x) dx.
% This is the direct trilayer analogue of the bilayer function plotted in
% Massatt-Carr-Luskin Figure 4.  It is not a single Hamiltonian coefficient.
%
% If includeBlochNormalization=true, the result is multiplied by
% sqrt(|Gamma_j^*||Gamma_k^*|), matching their tilde-h convention.
% Sublattice phases and chi(Q) are not included.
% No FFT is used.
if nargin < 8
    options = struct;
end
defaults = struct( ...
    'realSpaceRadialOrder',80, ...
    'realSpaceAngularOrder',128, ...
    'realSpaceCutoff',15*stack.aG, ...
    'configurationGridSize',[], ...
    'minimumConfigurationGridSize',9, ...
    'batchSize',128, ...
    'useGPU',false, ...
    'transformType','spectatorMode', ...
    'spectatorFraction',[0;0], ...
    'includeBlochNormalization',false);
options = mergePreserveEmpty(options,defaults,'configurationGridSize');

[j,k,l] = validatePair(pair);
alpha = alpha(:).';
beta = beta(:).';
assert(numel(alpha) == numel(beta) && ~isempty(alpha), ...
    'alpha and beta must have the same nonzero length.');
assert(all(ismember(alpha,1:2)) && all(ismember(beta,1:2)), ...
    'alpha and beta entries must be 1 or 2.');
assert(size(Q,1) == 2,'Q must be a 2-by-N array.');
spectatorMode = spectatorMode(:);
assert(numel(spectatorMode) == 2 && ...
    max(abs(spectatorMode-round(spectatorMode))) < 1e-12, ...
    'spectatorMode must be a two-component integer vector.');
spectatorMode = round(spectatorMode);

recommendedGridSize = max(options.minimumConfigurationGridSize, ...
    2*max(abs(spectatorMode))+1);
if isempty(options.configurationGridSize)
    configurationGridSize = recommendedGridSize;
else
    configurationGridSize = options.configurationGridSize;
    if configurationGridSize < recommendedGridSize
        warning('computeRelaxedInterlayerFourierHopping:aliasing', ...
            'configurationGridSize is below the recommended value %d.', ...
            recommendedGridSize);
    end
end

[x,xWeight] = polarQuadrature(options.realSpaceRadialOrder, ...
    options.realSpaceAngularOrder,options.realSpaceCutoff);
nX = size(x,2);
useGPU = options.useGPU;
hoppingMode = complex(zeros(nX,numel(alpha)));
switch lower(options.transformType)
    case {'spectatormode','mode'}
        [b,bWeight] = periodicCellTrapezoid( ...
            stack.A{l},configurationGridSize);
        nB = size(b,2);
        bDevice = toDevice(b,useGPU);
        bWeightDevice = toDevice(bWeight(:).',useGPU);
        dGlDevice = toDevice((stack.B{l}*spectatorMode).',useGPU);
        phaseB = exp(1i*(dGlDevice*bDevice)).*bWeightDevice;
        for channel = 1:numel(alpha)
            relaxedPosition = evaluateRelaxedInterlayerPosition( ...
                stack,relaxationFields,x,b,j,k,l, ...
                alpha(channel),beta(channel),useGPU);
            hopping = realSpaceInterlayerHopping( ...
                relaxedPosition,stack,j,k,alpha(channel),beta(channel));
            hopping = reshape(hopping,nX,nB);
            hoppingMode(:,channel) = toHost( ...
                hopping*phaseB.',useGPU);
        end
    case {'spectatorslice','slice'}
        fraction = options.spectatorFraction(:);
        assert(numel(fraction) == 2, ...
            'spectatorFraction must contain two components.');
        b = stack.A{l}*fraction;
        for channel = 1:numel(alpha)
            relaxedPosition = evaluateRelaxedInterlayerPosition( ...
                stack,relaxationFields,x,b,j,k,l, ...
                alpha(channel),beta(channel),useGPU);
            hopping = realSpaceInterlayerHopping( ...
                relaxedPosition,stack,j,k,alpha(channel),beta(channel));
            hoppingMode(:,channel) = toHost(hopping.',useGPU);
        end
    otherwise
        error('Unknown transformType: %s.',options.transformType);
end

xDevice = toDevice(x,useGPU);
xWeightDevice = toDevice(xWeight(:).',useGPU);
hoppingModeDevice = toDevice(hoppingMode,useGPU);
values = complex(zeros(size(Q,2),numel(alpha)));
for first = 1:options.batchSize:size(Q,2)
    last = min(first+options.batchSize-1,size(Q,2));
    range = first:last;
    QDevice = toDevice(Q(:,range).',useGPU);
    phaseX = exp(-1i*(QDevice*xDevice)).*xWeightDevice;
    values(range,:) = toHost( ...
        phaseX*hoppingModeDevice,useGPU)/(2*pi)^2;
end
if options.includeBlochNormalization
    values = values*sqrt(abs(det(stack.B{j}))*abs(det(stack.B{k})));
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

function options = mergePreserveEmpty(options,defaults,emptyName)
if isempty(options)
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
