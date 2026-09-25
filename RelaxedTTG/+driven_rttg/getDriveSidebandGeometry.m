function geometry = getDriveSidebandGeometry( ...
    q,opticalSources,acousticSources,hbar)
%GETDRIVESIDEBANDGEOMETRY Construct locked optical-acoustic sideband data.
if nargin < 4 || isempty(hbar)
    hbar = 1;
end
sources = driven_rttg.combineDriveSources( ...
    opticalSources,acousticSources);
numberOptical = numel(opticalSources);
numberAcoustic = numel(acousticSources);
numberSources = numel(sources);
assert(numberSources >= 1,'At least one drive source is required.');
assert(isequal(size(q),[2,1]),'q must be a 2-by-1 column vector.');

cutoff = zeros(1,numberSources);
waveVector = zeros(2,numberSources);
frequency = zeros(1,numberSources);
for source = 1:numberSources
    cutoff(source) = sources(source).cutoff;
    waveVector(:,source) = sources(source).waveVector;
    frequency(source) = sources(source).frequency;
end

axesValues = arrayfun(@(value) -value:value,cutoff, ...
    'UniformOutput',false);
grid = cell(1,numberSources);
[grid{:}] = ndgrid(axesValues{:});
numberSidebands = numel(grid{1});
sidebandIndex = zeros(numberSidebands,numberSources);
for source = 1:numberSources
    sidebandIndex(:,source) = grid{source}(:);
end

geometry.sidebandIndex = sidebandIndex;
geometry.opticalIndex = sidebandIndex(:,1:numberOptical);
geometry.acousticIndex = sidebandIndex(:,numberOptical+(1:numberAcoustic));
geometry.shiftedQ = q-waveVector*sidebandIndex.';
geometry.sidebandEnergy = hbar*(sidebandIndex*frequency.');
geometry.centralSideband = find(all(sidebandIndex == 0,2));
geometry.numberSidebands = numberSidebands;
geometry.numberOpticalSources = numberOptical;
geometry.numberAcousticSources = numberAcoustic;
geometry.waveVector = waveVector;
geometry.frequency = frequency;
geometry.cutoff = cutoff;
end
