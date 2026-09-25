function [figureHandle,data] = plotInterlayerMassattComparison( ...
    stack,relaxationFields,pair,options)
%PLOTINTERLAYERMASSATTCOMPARISON Compare unrelaxed/relaxed h and tilde h.
%
% The default Fourier transform is taken at a fixed spectator-layer slice,
% so it transforms exactly the function shown in the corresponding real-
% space panel.  This is the closest trilayer analogue of Figure 4 in
% Massatt, Carr, and Luskin.  The full trilayer Hamiltonian instead uses the
% spectatorMode coefficients; select that transform explicitly when needed.
if nargin < 4
    options = struct;
end
defaults = struct( ...
    'realSpace',struct, ...
    'fourier',struct, ...
    'qExtent',5, ...
    'qGridSize',101, ...
    'spectatorFraction',[0;0], ...
    'fourierComponent','magnitude', ...
    'channels',[1,1;1,2], ...
    'commonColorLimits',false);
options = mergeOptions(options,defaults);

options.realSpace = mergeOptions(options.realSpace,struct( ...
    'extent',5,'gridSize',161, ...
    'spectatorFraction',options.spectatorFraction,'useGPU',false));
options.fourier = mergeOptions(options.fourier,struct( ...
    'realSpaceRadialOrder',100, ...
    'realSpaceAngularOrder',180, ...
    'realSpaceCutoff',8, ...
    'configurationGridSize',31, ...
    'minimumConfigurationGridSize',9, ...
    'batchSize',128, ...
    'useGPU',false, ...
    'transformType','spectatorSlice', ...
    'spectatorFraction',options.spectatorFraction, ...
    'includeBlochNormalization',true));

assert(size(options.channels,2) == 2, ...
    'channels must have rows [alpha,beta].');
assert(size(options.channels,1) == 2, ...
    'This comparison layout expects exactly two channels.');
zeroFields = cell(1,3);
for layer = 1:3
    zeroFields{layer} = @(xk,xl) zeros(size(xk),'like',xk); %#ok<INUSD>
end

realData = cell(1,2);
realData{1} = sampleRelaxedInterlayerRealSpaceHopping( ...
    stack,zeroFields,pair,options.realSpace);
realData{2} = sampleRelaxedInterlayerRealSpaceHopping( ...
    stack,relaxationFields,pair,options.realSpace);

qAxis = linspace(-options.qExtent,options.qExtent,options.qGridSize);
[q1,q2] = ndgrid(qAxis,qAxis);
Q = [q1(:),q2(:)].';
fourierData = cell(2,2);
alphaList = options.channels(:,1).';
betaList = options.channels(:,2).';
for state = 1:2
    fields = zeroFields;
    if state == 2
        fields = relaxationFields;
    end
    values = computeRelaxedInterlayerFourierHopping( ...
        stack,fields,pair,alphaList,betaList,Q,[0;0],options.fourier);
    for channel = 1:2
        fourierData{state,channel} = reshape( ...
            values(:,channel),options.qGridSize,options.qGridSize);
    end
end

figureHandle = figure('Name','Interlayer hopping: Massatt comparison');
layout = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
realAxes = gobjects(2,2);
fourierAxes = gobjects(2,2);
for state = 1:2
    for channel = 1:2
        alpha = options.channels(channel,1);
        beta = options.channels(channel,2);
        realAxes(state,channel) = nexttile(layout, ...
            (state-1)*4+2*channel-1);
        imagesc(realData{state}.axis,realData{state}.axis, ...
            realData{state}.hopping{alpha,beta}.');
        axis xy equal tight
        xlabel('b_x (Angstrom)'); ylabel('b_y (Angstrom)');
        title(panelTitle(pair,alpha,beta,state,false),'Interpreter','latex');
        colorbar

        fourierAxes(state,channel) = nexttile(layout, ...
            (state-1)*4+2*channel);
        imagesc(qAxis,qAxis,selectComponent( ...
            fourierData{state,channel},options.fourierComponent).');
        axis xy equal tight
        xlabel('q_x (Angstrom^{-1})'); ylabel('q_y (Angstrom^{-1})');
        title(panelTitle(pair,alpha,beta,state,true),'Interpreter','latex');
        colorbar
    end
end

if options.commonColorLimits
    applyCommonLimits(realAxes);
    applyCommonLimits(fourierAxes);
end
sgtitle(layout,sprintf(['Interlayer pair (%d,%d), spectator slice ' ...
    '(%.3g,%.3g)'],pair(1),pair(2),options.spectatorFraction(1), ...
    options.spectatorFraction(2)));

data.realSpace = realData;
data.fourier.axis = qAxis;
data.fourier.values = fourierData;
data.fourier.transformType = options.fourier.transformType;
data.fourier.includeBlochNormalization = ...
    options.fourier.includeBlochNormalization;
data.options = options;
end

function titleText = panelTitle(pair,alpha,beta,state,isFourier)
labels = {'A','B'};
stateLabels = {'unrelaxed','relaxed'};
if isFourier
    symbol = '\tilde{h}';
else
    symbol = 'h';
end
titleText = sprintf('$%s_{%d%d}^{%s%s}$ %s',symbol,pair(1),pair(2), ...
    labels{alpha},labels{beta},stateLabels{state});
end

function value = selectComponent(z,component)
switch lower(component)
    case {'magnitude','abs'}
        value = abs(z);
    case 'real'
        value = real(z);
    case {'imaginary','imag'}
        value = imag(z);
    case 'phase'
        value = angle(z);
    otherwise
        error('Unknown fourierComponent: %s.',component);
end
end

function applyCommonLimits(axesHandles)
limits = cell2mat(arrayfun(@(ax) clim(ax),axesHandles(:), ...
    'UniformOutput',false));
common = [min(limits(:,1)),max(limits(:,2))];
for axisHandle = axesHandles(:).'
    clim(axisHandle,common);
end
end
