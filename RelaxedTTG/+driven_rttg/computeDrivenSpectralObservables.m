function observable = computeDrivenSpectralObservables( ...
    driveCache,E,eta,envelopeTime,options)
%COMPUTEDRIVENSPECTRALOBSERVABLES Central-sector driven spectral measures.
%
% This evaluates the central temporal-sector compression of
% (K_q(t)-z I)^(-1).  Its diagonal imaginary part gives the frozen-envelope
% version of the phase-averaged momentum spectral function in Eq. (42) of
% the supplied note.  No trace over temporal replicas is taken.
if nargin < 5
    options = struct;
end
defaults = struct( ...
    'spatialSelector',[], ...
    'lanczos',struct, ...
    'energyChunkSize',512, ...
    'computeOpticalAcousticContrast',false, ...
    'returnMatrixSpectralFunction',false, ...
    'verbose',true);
options = mergeOptions(options,defaults);

E = E(:).';
eta = eta(:);
envelopeTime = envelopeTime(:);
assert(~isempty(E) && all(isfinite(E)),'E must be finite and nonempty.');
assert(~isempty(eta) && all(isfinite(eta)) && all(eta > 0), ...
    'eta must contain positive finite broadenings.');
assert(~isempty(envelopeTime) && all(isfinite(envelopeTime)), ...
    'envelopeTime must be finite and nonempty.');
assert(isscalar(options.energyChunkSize) && ...
    options.energyChunkSize >= 1 && ...
    options.energyChunkSize == fix(options.energyChunkSize), ...
    'energyChunkSize must be a positive integer.');

if isempty(options.spatialSelector)
    assert(isfield(driveCache,'defaultSpatialSelector'), ...
        ['No default selector is stored in this drive cache. Supply ' ...
         'options.spatialSelector.']);
    spatialSelector = driveCache.defaultSpatialSelector;
else
    spatialSelector = options.spatialSelector;
end
selector = driven_rttg.getDrivenCentralSelector( ...
    driveCache,spatialSelector);
gram = full(selector'*selector);
assert(norm(gram-eye(size(gram)),'fro') < 1e-12, ...
    'The observable selector columns must be orthonormal.');
numberOrbitals = size(selector,2);

if options.computeOpticalAcousticContrast
    assert(~isempty(driveCache.opticalSources) && ...
        ~isempty(driveCache.acousticSources), ...
        'The optical-acoustic contrast requires both source families.');
end

numberTimes = numel(envelopeTime);
numberEta = numel(eta);
numberEnergy = numel(E);
observable.combined = allocateSpectralData( ...
    numberTimes,numberEta,numberEnergy,numberOrbitals, ...
    options.returnMatrixSpectralFunction);
if options.computeOpticalAcousticContrast
    observable.opticalOnly = allocateSpectralData( ...
        numberTimes,numberEta,numberEnergy,numberOrbitals,false);
    observable.acousticOnly = allocateSpectralData( ...
        numberTimes,numberEta,numberEnergy,numberOrbitals,false);
    observable.undriven = allocateSpectralData( ...
        numberTimes,numberEta,numberEnergy,numberOrbitals,false);
end

observable.diagnostics = cell(numberTimes,1);
for timeIndex = 1:numberTimes
    diagnostics = struct;
    [K,parts,info] = driven_rttg.evaluateDrivenSpaceTimeHamiltonian( ...
        driveCache,envelopeTime(timeIndex));
    combined = projectedSpectralData( ...
        K,selector,driveCache,E,eta,options);
    observable.combined = assignTime( ...
        observable.combined,combined,timeIndex);
    diagnostics.combined = combined.diagnostics;
    diagnostics.hamiltonian = info;

    if options.computeOpticalAcousticContrast
        opticalOnly = projectedSpectralData( ...
            parts.static+parts.optical,selector,driveCache,E,eta,options);
        acousticOnly = projectedSpectralData( ...
            parts.static+parts.acoustic,selector,driveCache,E,eta,options);
        undriven = projectedSpectralData( ...
            parts.static,selector,driveCache,E,eta,options);
        observable.opticalOnly = assignTime( ...
            observable.opticalOnly,opticalOnly,timeIndex);
        observable.acousticOnly = assignTime( ...
            observable.acousticOnly,acousticOnly,timeIndex);
        observable.undriven = assignTime( ...
            observable.undriven,undriven,timeIndex);
        diagnostics.opticalOnly = opticalOnly.diagnostics;
        diagnostics.acousticOnly = acousticOnly.diagnostics;
        diagnostics.undriven = undriven.diagnostics;
    end
    observable.diagnostics{timeIndex} = diagnostics;
    if options.verbose
        fprintf('Driven observables: completed envelope time %d of %d.\n', ...
            timeIndex,numberTimes);
    end
end

observable.combined = addLayerAverages(observable.combined);
if options.computeOpticalAcousticContrast
    observable.opticalOnly = addLayerAverages(observable.opticalOnly);
    observable.acousticOnly = addLayerAverages(observable.acousticOnly);
    observable.undriven = addLayerAverages(observable.undriven);
    observable.contrast.orbitalSpectralFunction = ...
        observable.combined.orbitalSpectralFunction ...
        -observable.opticalOnly.orbitalSpectralFunction ...
        -observable.acousticOnly.orbitalSpectralFunction ...
        +observable.undriven.orbitalSpectralFunction;
    observable.contrast.totalSpectralFunction = ...
        observable.combined.totalSpectralFunction ...
        -observable.opticalOnly.totalSpectralFunction ...
        -observable.acousticOnly.totalSpectralFunction ...
        +observable.undriven.totalSpectralFunction;
    observable.contrast.orbitalIntegratedMeasure = ...
        observable.combined.orbitalIntegratedMeasure ...
        -observable.opticalOnly.orbitalIntegratedMeasure ...
        -observable.acousticOnly.orbitalIntegratedMeasure ...
        +observable.undriven.orbitalIntegratedMeasure;
    observable.contrast.totalIntegratedMeasure = ...
        observable.combined.totalIntegratedMeasure ...
        -observable.opticalOnly.totalIntegratedMeasure ...
        -observable.acousticOnly.totalIntegratedMeasure ...
        +observable.undriven.totalIntegratedMeasure;
    if isfield(observable.combined,'layerSpectralFunction')
        observable.contrast.layerSpectralFunction = ...
            observable.combined.layerSpectralFunction ...
            -observable.opticalOnly.layerSpectralFunction ...
            -observable.acousticOnly.layerSpectralFunction ...
            +observable.undriven.layerSpectralFunction;
        observable.contrast.layerIntegratedMeasure = ...
            observable.combined.layerIntegratedMeasure ...
            -observable.opticalOnly.layerIntegratedMeasure ...
            -observable.acousticOnly.layerIntegratedMeasure ...
            +observable.undriven.layerIntegratedMeasure;
    end
end

observable.energy = E;
observable.eta = eta;
observable.envelopeTime = envelopeTime;
observable.selector = selector;
observable.orbitalLabels = defaultOrbitalLabels(numberOrbitals);
observable.sidebandIndex = driveCache.geometry.sidebandIndex;
observable.sidebandMomentumShift = driveCache.geometry.waveVector ...
    *driveCache.geometry.sidebandIndex.';
observable.sidebandEnergy = driveCache.geometry.sidebandEnergy;
observable.options = options;
observable.convention = ['Green=(K-zI)^(-1), spectralFunction=' ...
    'imag(Green)/pi, central temporal sector only'];
observable.interpretation = ['For pulsed drives this is a frozen-envelope ' ...
    'spectral diagnostic. A true time-resolved photoemission intensity ' ...
    'requires propagation, occupations, probe envelopes, and matrix ' ...
    'elements through a lesser Green function.'];
end

function data = projectedSpectralData(K,selector,driveCache,E,eta,options)
lanczos = rttg_common.blockLanczosHermitian( ...
    K,selector,options.lanczos);
nodes = lanczos.ritzValues(:);
orbitalWeight = abs(lanczos.greenCoupling).^2;
orbitalWeightError = max(abs(sum(orbitalWeight,2)-1));
assert(orbitalWeightError < 1e-10, ...
    'The projected pole weights do not have unit mass.');
numberOrbitals = size(selector,2);
numberEta = numel(eta);
numberEnergy = numel(E);

greenDiagonal = complex(zeros(numberEta,numberEnergy,numberOrbitals));
spectralFunction = zeros(numberEta,numberEnergy,numberOrbitals);
integratedMeasure = zeros(numberEta,numberEnergy,numberOrbitals);
if options.returnMatrixSpectralFunction
    matrixSpectralFunction = complex(zeros( ...
        numberEta,numberEnergy,numberOrbitals,numberOrbitals));
else
    matrixSpectralFunction = [];
end

for firstEnergy = 1:options.energyChunkSize:numberEnergy
    lastEnergy = min(firstEnergy+options.energyChunkSize-1,numberEnergy);
    range = firstEnergy:lastEnergy;
    energyBlock = E(range);
    for etaIndex = 1:numberEta
        inverseDenominator = 1./(nodes-(energyBlock+1i*eta(etaIndex)));
        diagonalBlock = orbitalWeight*inverseDenominator;
        greenDiagonal(etaIndex,range,:) = permute(diagonalBlock,[3,2,1]);
        spectralFunction(etaIndex,range,:) = ...
            permute(imag(diagonalBlock)/pi,[3,2,1]);
        cdfBlock = orbitalWeight*( ...
            0.5+atan((energyBlock-nodes)/eta(etaIndex))/pi);
        integratedMeasure(etaIndex,range,:) = permute(cdfBlock,[3,2,1]);

        if options.returnMatrixSpectralFunction
            coupling = lanczos.greenCoupling;
            for localEnergy = 1:numel(range)
                greenMatrix = (coupling ...
                    .*inverseDenominator(:,localEnergy).')*coupling';
                densityMatrix = (greenMatrix-greenMatrix')/(2i*pi);
                matrixSpectralFunction( ...
                    etaIndex,range(localEnergy),:,:) = reshape( ...
                    densityMatrix,1,1,numberOrbitals,numberOrbitals);
            end
        end
    end
end

boundarySideband = false(driveCache.geometry.numberSidebands,1);
for source = 1:numel(driveCache.geometry.cutoff)
    cutoff = driveCache.geometry.cutoff(source);
    if cutoff > 0
        boundarySideband = boundarySideband ...
            |abs(driveCache.geometry.sidebandIndex(:,source)) == cutoff;
    end
end
boundaryRows = sidebandRows( ...
    find(boundarySideband),driveCache.baseDimension);
centralRows = driveCache.centralSpatialRange;
boundaryAmplitude = lanczos.basis(boundaryRows,:)*lanczos.ritzVectors;
centralAmplitude = lanczos.basis(centralRows,:)*lanczos.ritzVectors;

data.greenDiagonal = greenDiagonal;
data.orbitalSpectralFunction = spectralFunction;
data.orbitalIntegratedMeasure = integratedMeasure;
data.matrixSpectralFunction = matrixSpectralFunction;
data.padeNodes = nodes;
data.padeOrbitalWeights = orbitalWeight;
data.boundaryWeightByPole = sum(abs(boundaryAmplitude).^2,1).';
data.centralWeightByPole = sum(abs(centralAmplitude).^2,1).';
scalarPoleWeight = mean(orbitalWeight,1).';
data.diagnostics.orthogonalityError = lanczos.orthogonalityError;
data.diagnostics.blockTridiagonalDefect = ...
    lanczos.blockTridiagonalDefect;
data.diagnostics.krylovDimension = lanczos.krylovDimension;
data.diagnostics.numberBlockSteps = lanczos.numberBlockSteps;
data.diagnostics.orbitalWeightError = orbitalWeightError;
data.diagnostics.maximumBoundaryWeight = ...
    max(data.boundaryWeightByPole);
data.diagnostics.centralMeasureWeightedBoundaryWeight = sum( ...
    scalarPoleWeight.*data.boundaryWeightByPole);
end

function rows = sidebandRows(sideband,baseDimension)
rows = reshape((sideband(:)-1)*baseDimension+(1:baseDimension),[],1);
end

function data = allocateSpectralData( ...
    numberTimes,numberEta,numberEnergy,numberOrbitals,includeMatrix)
data.greenDiagonal = complex(zeros( ...
    numberTimes,numberEta,numberEnergy,numberOrbitals));
data.orbitalSpectralFunction = zeros( ...
    numberTimes,numberEta,numberEnergy,numberOrbitals);
data.orbitalIntegratedMeasure = zeros( ...
    numberTimes,numberEta,numberEnergy,numberOrbitals);
if includeMatrix
    data.matrixSpectralFunction = complex(zeros( ...
        numberTimes,numberEta,numberEnergy,numberOrbitals,numberOrbitals));
else
    data.matrixSpectralFunction = [];
end
data.padeNodes = cell(numberTimes,1);
data.padeOrbitalWeights = cell(numberTimes,1);
data.boundaryWeightByPole = cell(numberTimes,1);
data.centralWeightByPole = cell(numberTimes,1);
end

function destination = assignTime(destination,source,timeIndex)
numberEta = size(source.greenDiagonal,1);
numberEnergy = size(source.greenDiagonal,2);
numberOrbitals = size(source.greenDiagonal,3);
destination.greenDiagonal(timeIndex,:,:,:) = reshape( ...
    source.greenDiagonal,1,numberEta,numberEnergy,numberOrbitals);
destination.orbitalSpectralFunction(timeIndex,:,:,:) = ...
    reshape(source.orbitalSpectralFunction, ...
    1,numberEta,numberEnergy,numberOrbitals);
destination.orbitalIntegratedMeasure(timeIndex,:,:,:) = ...
    reshape(source.orbitalIntegratedMeasure, ...
    1,numberEta,numberEnergy,numberOrbitals);
if ~isempty(destination.matrixSpectralFunction)
    destination.matrixSpectralFunction(timeIndex,:,:,:,:) = ...
        reshape(source.matrixSpectralFunction, ...
        1,numberEta,numberEnergy,numberOrbitals,numberOrbitals);
end
destination.padeNodes{timeIndex} = source.padeNodes;
destination.padeOrbitalWeights{timeIndex} = source.padeOrbitalWeights;
destination.boundaryWeightByPole{timeIndex} = source.boundaryWeightByPole;
destination.centralWeightByPole{timeIndex} = source.centralWeightByPole;
end

function data = addLayerAverages(data)
data.totalSpectralFunction = mean(data.orbitalSpectralFunction,4);
data.totalIntegratedMeasure = mean(data.orbitalIntegratedMeasure,4);
if size(data.orbitalSpectralFunction,4) == 6
    data.layerSpectralFunction = zeros( ...
        size(data.orbitalSpectralFunction,1), ...
        size(data.orbitalSpectralFunction,2), ...
        size(data.orbitalSpectralFunction,3),3);
    data.layerIntegratedMeasure = zeros(size(data.layerSpectralFunction));
    for layer = 1:3
        columns = 2*layer+(-1:0);
        data.layerSpectralFunction(:,:,:,layer) = mean( ...
            data.orbitalSpectralFunction(:,:,:,columns),4);
        data.layerIntegratedMeasure(:,:,:,layer) = mean( ...
            data.orbitalIntegratedMeasure(:,:,:,columns),4);
    end
end
end

function labels = defaultOrbitalLabels(numberOrbitals)
if numberOrbitals == 6
    labels = {'1A','1B','2A','2B','3A','3B'};
else
    labels = compose('orbital %d',1:numberOrbitals);
end
end
