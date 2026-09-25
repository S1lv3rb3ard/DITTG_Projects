function [figureHandle,axesHandle,surfaceHandle] = ...
    plotMomentumLocalIDOSLinecut( ...
        linecut,E,momentumLocalIDOS,eta,options)
%PLOTMOMENTUMLOCALIDOSLINECUT Plot N_eta(q,E) over a symmetry line cut.
%
% momentumLocalIDOS must be numberQ-by-numberEnergy or
% numberQ-by-numberEnergy-by-numberEta.  options.etaIndex selects the
% smoothing width when a three-dimensional array is supplied.
if nargin < 5
    options = struct;
end
defaults = struct( ...
    'etaIndex',1, ...
    'figureVisible','on', ...
    'title','Momentum-local integrated spectral measure', ...
    'colorLimits',[], ...
    'colormap','turbo', ...
    'view',[45,32]);
options = mergeOptions(options,defaults);

E = E(:).';
eta = eta(:);
numberQ = numel(linecut.kPath);
assert(size(linecut.q,2) == numberQ, ...
    'linecut.q and linecut.kPath must have the same number of points.');
assert(numel(E) >= 2 && numberQ >= 2, ...
    'At least two energies and two path points are required.');
assert(options.etaIndex >= 1 && options.etaIndex <= numel(eta) && ...
    options.etaIndex == fix(options.etaIndex), ...
    'options.etaIndex must select an entry of eta.');

if ismatrix(momentumLocalIDOS)
    surfaceData = momentumLocalIDOS;
    assert(numel(eta) == 1 || options.etaIndex == 1, ...
        'A two-dimensional IDOS array contains only one eta slice.');
else
    assert(size(momentumLocalIDOS,3) == numel(eta), ...
        'The third IDOS dimension must agree with eta.');
    surfaceData = momentumLocalIDOS(:,:,options.etaIndex);
end
assert(isequal(size(surfaceData),[numberQ,numel(E)]), ...
    'IDOS must be numberQ-by-numberEnergy-by-numberEta.');

[pathGrid,energyGrid] = meshgrid(linecut.kPath,E);
figureHandle = figure('Visible',options.figureVisible,'Color','w');
axesHandle = axes(figureHandle);
surfaceHandle = surf(axesHandle,pathGrid,energyGrid, ...
    real(surfaceData.'),real(surfaceData.'), ...
    'EdgeColor','none','FaceColor','interp');
view(axesHandle,options.view);
xlim(axesHandle,[linecut.kPath(1),linecut.kPath(end)]);
xticks(axesHandle,linecut.tickLocs);
xticklabels(axesHandle,linecut.labels);
set(axesHandle,'TickLabelInterpreter','latex');
xlabel(axesHandle,'high-symmetry path');
ylabel(axesHandle,'energy (eV)');
zlabel(axesHandle,'N_\eta(q,E)','Interpreter','tex');
title(axesHandle,sprintf('%s, \\eta = %.4g eV', ...
    options.title,eta(options.etaIndex)),'Interpreter','tex');
colormap(axesHandle,options.colormap);
colorbar(axesHandle);
grid(axesHandle,'on');
if ~isempty(options.colorLimits)
    assert(numel(options.colorLimits) == 2 && ...
        options.colorLimits(2) > options.colorLimits(1), ...
        'options.colorLimits must be [minimum maximum].');
    caxis(axesHandle,options.colorLimits);
end
end
