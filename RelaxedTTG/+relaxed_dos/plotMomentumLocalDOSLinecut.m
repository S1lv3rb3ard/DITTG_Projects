function [figureHandle,axesHandle,imageHandle] = ...
    plotMomentumLocalDOSLinecut(linecut,E,LDoS,options)
%PLOTMOMENTUMLDOSLINECUT Plot momentum LDoS versus path distance and energy.
%
% LDoS must be numberQ-by-numberEnergy, as returned in result.LDoS by
% relaxed_dos.computeMomentumLocalDOSLinecut.
if nargin < 4
    options = struct;
end
defaults = struct( ...
    'figureVisible','on', ...
    'title','Relaxed momentum LDoS', ...
    'colorLimits',[], ...
    'colormap','turbo');
options = mergeOptions(options,defaults);

E = E(:).';
assert(size(linecut.q,2) == numel(linecut.kPath), ...
    'linecut.q and linecut.kPath must have the same number of points.');
assert(numel(linecut.kPath) >= 2 && numel(E) >= 2, ...
    'At least two path points and two energies are required for plotting.');
assert(isequal(size(LDoS),[numel(linecut.kPath),numel(E)]), ...
    'LDoS must be numberQ-by-numberEnergy.');

figureHandle = figure('Visible',options.figureVisible,'Color','w');
axesHandle = axes(figureHandle);
[pathGrid,energyGrid] = meshgrid(linecut.kPath,E);
imageHandle = surface(axesHandle,pathGrid,energyGrid, ...
    zeros(size(energyGrid)),real(LDoS.'), ...
    'EdgeColor','none','FaceColor','interp');
view(axesHandle,2);
set(axesHandle,'YDir','normal');
xlim(axesHandle,[linecut.kPath(1),linecut.kPath(end)]);
xticks(axesHandle,linecut.tickLocs);
xticklabels(axesHandle,linecut.labels);
set(axesHandle,'TickLabelInterpreter','latex','Layer','top');
xlabel(axesHandle,'high-symmetry path');
ylabel(axesHandle,'energy (eV)');
title(axesHandle,options.title,'Interpreter','none');
colormap(axesHandle,options.colormap);
colorbar(axesHandle);
if ~isempty(options.colorLimits)
    assert(numel(options.colorLimits) == 2 && ...
        options.colorLimits(2) > options.colorLimits(1), ...
        'options.colorLimits must be [minimum maximum].');
    caxis(axesHandle,options.colorLimits);
end
end
