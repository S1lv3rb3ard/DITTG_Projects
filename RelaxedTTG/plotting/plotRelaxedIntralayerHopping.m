function [figureHandle,data] = plotRelaxedIntralayerHopping( ...
    stack,shellsReference,relaxationFields,layer,options)
%PLOTRELAXEDINTRALAYERHOPPING Plot relaxed discrete bond hoppings for a layer.
%
% A fixed configuration slice is selected by fractional coordinates in the
% two complementary layer cells. Each panel shows the relaxed physical bond
% vectors, colored by hopping strength, for AA, AB, BA, and BB channels.
if nargin < 5
    options = struct;
end
defaults = struct( ...
    'configurationFractionK',[0;0], ...
    'configurationFractionL',[0;0], ...
    'shellReferenceLayer',2, ...
    'markerSize',60, ...
    'figureHandle',[]);
options = mergeOptions(options,defaults);

assert(ismember(layer,1:3),'layer must be 1, 2, or 3.');
assert(iscell(relaxationFields) && numel(relaxationFields) == 3, ...
    'relaxationFields must contain three function handles.');
fractionK = options.configurationFractionK(:);
fractionL = options.configurationFractionL(:);
assert(numel(fractionK) == 2 && numel(fractionL) == 2, ...
    'Configuration fractions must contain two components.');

otherLayers = setdiff(1:3,layer,'stable');
k = otherLayers(1);
l = otherLayers(2);
bK = stack.A{k}*fractionK;
bL = stack.A{l}*fractionL;
t = stack.tau(:,layer);
assert(ismember(options.shellReferenceLayer,1:3), ...
    'shellReferenceLayer must be 1, 2, or 3.');
shellMap = stack.A{layer}/stack.A{options.shellReferenceLayer};
RA = shellMap*shellsReference.matA;
RB = shellMap*shellsReference.matB;
u = relaxationFields{layer};

data = cell(2,2);
data{1,1} = makeChannel(RA, ...
    RA+u(RA+bK,RA+bL)-u(bK,bL),shellsReference.intraAA,'AA');
data{2,2} = makeChannel(RA, ...
    RA+u(RA+bK+t,RA+bL+t)-u(bK+t,bL+t), ...
    shellsReference.intraAA,'BB');
data{1,2} = makeChannel(RB, ...
    RB+u(RB+bK,RB+bL)-u(bK+t,bL+t), ...
    shellsReference.intraAB,'AB');
data{2,1} = makeChannel(RB, ...
    RB+u(RB+bK+t,RB+bL+t)-u(bK,bL), ...
    shellsReference.intraAB,'BA');

if isempty(options.figureHandle)
    figureHandle = figure('Name',sprintf('Relaxed intralayer %d',layer));
else
    figureHandle = options.figureHandle;
    figure(figureHandle);
    clf(figureHandle);
end
layout = tiledlayout(figureHandle,2,2,'TileSpacing','compact');
labels = {'AA','AB';'BA','BB'};
for alpha = 1:2
    for beta = 1:2
        nexttile(layout);
        channel = data{alpha,beta};
        scatter(channel.relaxedVector(1,:),channel.relaxedVector(2,:), ...
            options.markerSize,channel.hopping,'filled');
        hold on
        plot(0,0,'kx','LineWidth',1.5,'MarkerSize',8);
        hold off
        axis equal
        grid on
        colorbar
        xlabel('R_x');
        ylabel('R_y');
        title(sprintf('%s hopping',labels{alpha,beta}));
    end
end
sgtitle(layout,sprintf('Layer %d relaxed intralayer bonds',layer));
end

function channel = makeChannel(unrelaxedVector,relaxedVector,interpolant,label)
channel.label = label;
channel.unrelaxedVector = unrelaxedVector;
channel.relaxedVector = relaxedVector;
channel.relaxedLength = vecnorm(relaxedVector,2,1);
channel.hopping = interpolant(channel.relaxedLength).';
channel.hopping = channel.hopping(:).';
end
