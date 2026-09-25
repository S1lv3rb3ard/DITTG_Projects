function report = checkDriveIncommensurability( ...
    opticalSources,acousticSources,options)
%CHECKDRIVEINCOMMENSURABILITY Search for low-order integer relations.
%
% This is a finite diagnostic, not a proof of rational independence. It
% searches |coefficient| <= relationSearchOrder for a simultaneous near-zero
% spatial and temporal combination of all optical and acoustic carriers.
if nargin < 3
    options = struct;
end
defaults = struct( ...
    'relationSearchOrder',4, ...
    'relationTolerance',1e-10, ...
    'maximumRelationCandidates',1e6);
options = mergeOptions(options,defaults);
assert(options.relationSearchOrder >= 1 && ...
    options.relationSearchOrder == fix(options.relationSearchOrder), ...
    'relationSearchOrder must be a positive integer.');
assert(options.relationTolerance > 0, ...
    'relationTolerance must be positive.');
sources = driven_rttg.combineDriveSources( ...
    opticalSources,acousticSources);
numberSources = numel(sources);
assert(numberSources >= 1,'At least one drive source is required.');

waveVector = zeros(2,numberSources);
frequency = zeros(1,numberSources);
for source = 1:numberSources
    waveVector(:,source) = sources(source).waveVector;
    frequency(source) = sources(source).frequency;
end

order = options.relationSearchOrder;
numberCandidates = (2*order+1)^numberSources-1;
assert(numberCandidates <= options.maximumRelationCandidates, ...
    ['Integer-relation search has %g candidates. Reduce ' ...
     'relationSearchOrder or increase maximumRelationCandidates.'], ...
    numberCandidates);

axesValues = repmat({-order:order},1,numberSources);
grid = cell(1,numberSources);
[grid{:}] = ndgrid(axesValues{:});
coefficient = zeros(numberSources,numberCandidates+1);
for source = 1:numberSources
    coefficient(source,:) = grid{source}(:).';
end
coefficient(:,all(coefficient == 0,1)) = [];

spatialResidual = waveVector*coefficient;
frequencyResidual = frequency*coefficient;
absoluteCoefficient = abs(coefficient);
spatialDenominator = vecnorm(waveVector,2,1)*absoluteCoefficient;
frequencyDenominator = abs(frequency)*absoluteCoefficient;
relativeSpatial = vecnorm(spatialResidual,2,1) ...
    ./max(spatialDenominator,eps);
relativeFrequency = abs(frequencyResidual) ...
    ./max(frequencyDenominator,eps);
jointResidual = max(relativeSpatial,relativeFrequency);

[minimumResidual,location] = min(jointResidual);
[minimumFrequencyResidual,frequencyLocation] = min(relativeFrequency);
report.minimumRelativeResidual = minimumResidual;
report.closestIntegerRelation = coefficient(:,location).';
report.spatialResidual = spatialResidual(:,location);
report.frequencyResidual = frequencyResidual(location);
report.detectedRelation = minimumResidual < options.relationTolerance;
report.minimumRelativeFrequencyResidual = minimumFrequencyResidual;
report.closestFrequencyRelation = coefficient(:,frequencyLocation).';
report.detectedFrequencyRelation = ...
    minimumFrequencyResidual < options.relationTolerance;
report.searchOrder = order;
report.numberCandidates = size(coefficient,2);
report.note = ['A finite integer search can detect low-order resonances ' ...
    'but cannot prove exact incommensurability.'];
end
