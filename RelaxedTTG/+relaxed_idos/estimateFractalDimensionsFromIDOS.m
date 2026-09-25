function dimension = estimateFractalDimensionsFromIDOS( ...
    E,IDOS,scales,options)
%ESTIMATEFRACTALDIMENSIONSFROMIDOS Estimate measure dimensions from IDoS.
%
% Local dimensions use
%   N(E0+r)-N(E0-r) ~ r^d.
% Generalized dimensions use box masses obtained from IDoS increments.
% IDOS may have one row reused at every scale, or one row per scale.  In
% the latter case options.smoothingWidths records the Green broadening used
% for each row and permits exclusion of scales too close to the smoothing.
if nargin < 4
    options = struct;
end
defaults = struct( ...
    'evaluationEnergies',[], ...
    'qOrders',[0,1,2], ...
    'window',[], ...
    'fitScaleIndices',[], ...
    'smoothingWidths',[], ...
    'minimumScaleToSmoothingRatio',2, ...
    'minimumResolvedScale',0, ...
    'minimumMass',1e-14, ...
    'occupiedMassTolerance',1e-12, ...
    'enforceMonotonicity',true);
options = mergeOptions(options,defaults);

assert(isscalar(options.minimumScaleToSmoothingRatio) && ...
    options.minimumScaleToSmoothingRatio >= 0, ...
    'minimumScaleToSmoothingRatio must be a nonnegative scalar.');
assert(isscalar(options.minimumResolvedScale) && ...
    options.minimumResolvedScale >= 0, ...
    'minimumResolvedScale must be a nonnegative scalar.');
assert(isvector(options.qOrders) && all(isfinite(options.qOrders)), ...
    'qOrders must be a finite vector.');

E = E(:).';
assert(numel(E) >= 3 && all(isfinite(E)) && all(diff(E) > 0), ...
    'E must be a strictly increasing finite vector with at least 3 points.');
scales = scales(:);
assert(numel(scales) >= 3 && all(isfinite(scales)) && all(scales > 0), ...
    'At least three finite positive scales are required.');

if isvector(IDOS)
    IDOS = IDOS(:).';
elseif size(IDOS,2) ~= numel(E) && size(IDOS,1) == numel(E)
    IDOS = IDOS.';
end
assert(size(IDOS,2) == numel(E) && ...
    any(size(IDOS,1) == [1,numel(scales)]), ...
    'IDOS must have one row or one row per scale.');
assert(all(isfinite(IDOS),'all'),'IDOS must contain only finite values.');
if size(IDOS,1) == 1
    IDOS = repmat(IDOS,numel(scales),1);
end

if options.enforceMonotonicity
    IDOS = cummax(IDOS,2);
else
    assert(all(diff(IDOS,1,2) >= -1e-12,'all'), ...
        'Each IDoS row must be nondecreasing.');
end

if isempty(options.window)
    window = [E(1),E(end)];
else
    window = options.window(:).';
end
assert(numel(window) == 2 && window(1) < window(2) && ...
    window(1) >= E(1) && window(2) <= E(end), ...
    'options.window must lie inside the sampled energy interval.');

if isempty(options.evaluationEnergies)
    evaluationEnergies = mean(window);
else
    evaluationEnergies = options.evaluationEnergies(:).';
end
assert(all(evaluationEnergies > window(1) & ...
    evaluationEnergies < window(2)), ...
    'Local-dimension evaluation energies must lie inside the window.');

if isempty(options.smoothingWidths)
    smoothingWidths = zeros(size(scales));
else
    smoothingWidths = options.smoothingWidths(:);
    assert(numel(smoothingWidths) == numel(scales) && ...
        all(smoothingWidths >= 0), ...
        'smoothingWidths must contain one nonnegative value per scale.');
end
if isempty(options.fitScaleIndices)
    fitScaleIndices = (1:numel(scales)).';
else
    fitScaleIndices = options.fitScaleIndices(:);
    assert(all(fitScaleIndices >= 1 & ...
        fitScaleIndices <= numel(scales) & ...
        fitScaleIndices == fix(fitScaleIndices)), ...
        'fitScaleIndices contains an invalid scale index.');
end

numberLocal = numel(evaluationEnergies);
localMass = NaN(numberLocal,numel(scales));
localDimension = NaN(1,numberLocal);
localFitR2 = NaN(1,numberLocal);
localEffectiveSlope = NaN(numberLocal,numel(scales)-1);
localLowerEstimate = NaN(1,numberLocal);
localUpperEstimate = NaN(1,numberLocal);
for energyIndex = 1:numberLocal
    center = evaluationEnergies(energyIndex);
    for scaleIndex = 1:numel(scales)
        radius = scales(scaleIndex);
        if center-radius < E(1) || center+radius > E(end)
            continue
        end
        row = IDOS(scaleIndex,:);
        upper = interp1(E,row,center+radius,'linear');
        lower = interp1(E,row,center-radius,'linear');
        localMass(energyIndex,scaleIndex) = max(upper-lower,0);
    end
    valid = false(numel(scales),1);
    valid(fitScaleIndices) = true;
    valid = valid & isfinite(localMass(energyIndex,:)).' ...
        & localMass(energyIndex,:).' > options.minimumMass ...
        & scales >= options.minimumScaleToSmoothingRatio*smoothingWidths ...
        & scales >= options.minimumResolvedScale;
    [localDimension(energyIndex),~,localFitR2(energyIndex)] = ...
        logarithmicFit(scales(valid),localMass(energyIndex,valid).');
    localEffectiveSlope(energyIndex,:) = diff(log(localMass(energyIndex,:))) ...
        ./diff(log(scales.'));
    pairValid = valid(1:end-1) & valid(2:end) ...
        & isfinite(localEffectiveSlope(energyIndex,:)).';
    if any(pairValid)
        localLowerEstimate(energyIndex) = min( ...
            localEffectiveSlope(energyIndex,pairValid));
        localUpperEstimate(energyIndex) = max( ...
            localEffectiveSlope(energyIndex,pairValid));
    end
end

qOrders = options.qOrders(:).';
numberQ = numel(qOrders);
partitionStatistic = NaN(numberQ,numel(scales));
actualBoxWidth = NaN(numel(scales),1);
numberOccupied = NaN(numel(scales),1);
for scaleIndex = 1:numel(scales)
    numberBoxes = max(1,round(diff(window)/scales(scaleIndex)));
    edges = linspace(window(1),window(2),numberBoxes+1);
    actualBoxWidth(scaleIndex) = edges(2)-edges(1);
    cdf = interp1(E,IDOS(scaleIndex,:),edges,'linear');
    probability = max(diff(cdf),0);
    totalMass = sum(probability);
    if totalMass <= options.minimumMass
        continue
    end
    probability = probability/totalMass;
    occupied = probability > ...
        options.occupiedMassTolerance*max(probability);
    positiveProbability = probability(occupied);
    numberOccupied(scaleIndex) = nnz(occupied);
    for qIndex = 1:numberQ
        q = qOrders(qIndex);
        if q == 0
            partitionStatistic(qIndex,scaleIndex) = ...
                log(numberOccupied(scaleIndex));
        elseif q == 1
            partitionStatistic(qIndex,scaleIndex) = ...
                sum(positiveProbability.*log(positiveProbability));
        else
            partitionStatistic(qIndex,scaleIndex) = ...
                log(sum(positiveProbability.^q));
        end
    end
end

generalizedDimension = NaN(size(qOrders));
generalizedFitR2 = NaN(size(qOrders));
generalizedEffectiveDimension = NaN(numberQ,numel(scales)-1);
for qIndex = 1:numberQ
    valid = false(numel(scales),1);
    valid(fitScaleIndices) = true;
    valid = valid & isfinite(partitionStatistic(qIndex,:)).' ...
        & actualBoxWidth >= ...
        options.minimumScaleToSmoothingRatio*smoothingWidths ...
        & actualBoxWidth >= options.minimumResolvedScale;
    [slope,~,generalizedFitR2(qIndex)] = logarithmicFit( ...
        actualBoxWidth(valid),partitionStatistic(qIndex,valid).',false);
    q = qOrders(qIndex);
    if q == 0
        generalizedDimension(qIndex) = -slope;
    elseif q == 1
        generalizedDimension(qIndex) = slope;
    else
        generalizedDimension(qIndex) = slope/(q-1);
    end
    effectiveSlope = diff(partitionStatistic(qIndex,:)) ...
        ./diff(log(actualBoxWidth.'));
    if q == 0
        generalizedEffectiveDimension(qIndex,:) = -effectiveSlope;
    elseif q == 1
        generalizedEffectiveDimension(qIndex,:) = effectiveSlope;
    else
        generalizedEffectiveDimension(qIndex,:) = effectiveSlope/(q-1);
    end
end

dimension.energy = E;
dimension.scales = scales;
dimension.smoothingWidths = smoothingWidths;
dimension.minimumResolvedScale = options.minimumResolvedScale;
dimension.window = window;
dimension.local.energy = evaluationEnergies;
dimension.local.mass = localMass;
dimension.local.dimension = localDimension;
dimension.local.fitR2 = localFitR2;
dimension.local.effectiveSlope = localEffectiveSlope;
dimension.local.lowerEstimate = localLowerEstimate;
dimension.local.upperEstimate = localUpperEstimate;
dimension.generalized.q = qOrders;
dimension.generalized.dimension = generalizedDimension;
dimension.generalized.fitR2 = generalizedFitR2;
dimension.generalized.effectiveDimension = generalizedEffectiveDimension;
dimension.generalized.partitionStatistic = partitionStatistic;
dimension.generalized.actualBoxWidth = actualBoxWidth;
dimension.generalized.numberOccupiedBoxes = numberOccupied;
dimension.warning = ['These are finite-scale regression estimates. ' ...
    'A fractal dimension requires a stable joint limit in truncation, ' ...
    'Lanczos depth, q quadrature, smoothing width, and fitting scale.'];
end

function [slope,intercept,R2] = logarithmicFit(x,y,takeLogY)
if nargin < 3
    takeLogY = true;
end
valid = isfinite(x) & isfinite(y) & x > 0;
if takeLogY
    valid = valid & y > 0;
end
x = x(valid);
y = y(valid);
if numel(x) < 3
    slope = NaN;
    intercept = NaN;
    R2 = NaN;
    return
end
logX = log(x);
if takeLogY
    fitY = log(y);
else
    fitY = y;
end
coefficient = polyfit(logX,fitY,1);
slope = coefficient(1);
intercept = coefficient(2);
residual = fitY-polyval(coefficient,logX);
R2 = 1-sum(residual.^2)/max(sum((fitY-mean(fitY)).^2),eps);
end
