function [figures,data] = plotRelaxedInterlayerHopping( ...
    stack,relaxationFields,pair,options)
%PLOTRELAXEDINTERLAYERHOPPING Plot real and Fourier-space hopping channels.
if nargin < 4
    options = struct;
end
defaults = struct( ...
    'realSpace',struct, ...
    'fourier',struct, ...
    'qExtent',2*pi, ...
    'qGridSize',61, ...
    'spectatorMode',[0;0], ...
    'fourierComponent','magnitude');
options = mergeOptions(options,defaults);
assert(options.qGridSize >= 2 && options.qGridSize == fix(options.qGridSize), ...
    'qGridSize must be an integer at least two.');

realData = sampleRelaxedInterlayerRealSpaceHopping( ...
    stack,relaxationFields,pair,options.realSpace);
j = realData.layers(1);
k = realData.layers(2);
l = realData.spectatorLayer;

figures.realSpace = figure('Name', ...
    sprintf('Relaxed interlayer real space %d-%d',j,k));
layout = tiledlayout(figures.realSpace,2,2,'TileSpacing','compact');
for alpha = 1:2
    for beta = 1:2
        nexttile(layout);
        imagesc(realData.axis,realData.axis, ...
            realData.hopping{alpha,beta}.');
        axis xy equal tight
        colorbar
        xlabel('x_1');
        ylabel('x_2');
        title(channelTitle(j,k,alpha,beta));
    end
end
sgtitle(layout,sprintf('Relaxed real-space hopping, b_%d slice',l));

qAxis = linspace(-options.qExtent,options.qExtent,options.qGridSize);
[q1,q2] = ndgrid(qAxis,qAxis);
Q = [q1(:),q2(:)].';
fourierData = cell(2,2);

figures.fourier = figure('Name', ...
    sprintf('Relaxed interlayer Fourier space %d-%d',j,k));
layout = tiledlayout(figures.fourier,2,2,'TileSpacing','compact');
for alpha = 1:2
    for beta = 1:2
        values = computeRelaxedInterlayerFourierHopping( ...
            stack,relaxationFields,pair,alpha,beta,Q, ...
            options.spectatorMode,options.fourier);
        fourierData{alpha,beta} = reshape( ...
            values,options.qGridSize,options.qGridSize);
        nexttile(layout);
        imagesc(qAxis,qAxis,selectComponent( ...
            fourierData{alpha,beta},options.fourierComponent).');
        axis xy equal tight
        colorbar
        xlabel('Q_1');
        ylabel('Q_2');
        title(channelTitle(j,k,alpha,beta));
    end
end
mode = options.spectatorMode(:);
if isfield(options.fourier,'transformType') && ...
        ismember(lower(options.fourier.transformType),{'spectatorslice','slice'})
    if isfield(options.fourier,'spectatorFraction')
        fraction = options.fourier.spectatorFraction(:);
    else
        fraction = [0;0];
    end
    sgtitle(layout,sprintf(['Relaxed Fourier hopping of spectator-layer ' ...
        '%d slice (%.3g,%.3g), %s'],l,fraction(1),fraction(2), ...
        options.fourierComponent));
else
    sgtitle(layout,sprintf(['Relaxed Fourier hopping, spectator layer %d, ' ...
        'mode (%d,%d), %s'],l,mode(1),mode(2),options.fourierComponent));
end

data.realSpace = realData;
data.fourier.axis = qAxis;
data.fourier.q1 = q1;
data.fourier.q2 = q2;
data.fourier.values = fourierData;
data.fourier.spectatorMode = mode;
data.fourier.component = options.fourierComponent;
if isfield(options.fourier,'transformType')
    data.fourier.transformType = options.fourier.transformType;
else
    data.fourier.transformType = 'spectatorMode';
end
end

function titleText = channelTitle(j,k,alpha,beta)
labels = {'A','B'};
titleText = sprintf('h_{%d%s}^{%d%s}', ...
    k,labels{beta},j,labels{alpha});
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
