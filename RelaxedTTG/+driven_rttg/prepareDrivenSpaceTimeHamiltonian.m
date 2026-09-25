function cache = prepareDrivenSpaceTimeHamiltonian( ...
    Hamiltonian,q,baseDimension,opticalSources,acousticSources,options)
%PREPAREDRIVENSPACETIMEHAMILTONIAN Prepare the locked extended Hamiltonian.
%
% The diagonal sideband indexed by (n,m) contains
%
%   H(q-n*K-m*Q) + hbar*(n*omega+m*Omega) I.
%
% Optical pulse envelopes are evaluated later by
% evaluateDrivenSpaceTimeHamiltonian, making the extended operator K_q(t)
% explicitly time dependent while retaining the carrier sideband indices.
if nargin < 6
    options = struct;
end
defaults = struct( ...
    'hbar',1, ...
    'maximumExtendedDimension',2e5, ...
    'maximumOpticalDutyCycle',0.10, ...
    'requirePulsedOpticalSources',true, ...
    'requireNoDetectedRelations',true, ...
    'requireNoDetectedFrequencyRelations',true, ...
    'relationSearchOrder',4, ...
    'relationTolerance',1e-10, ...
    'maximumRelationCandidates',1e6, ...
    'hermiticityTolerance',1e-12, ...
    'verbose',true);
options = mergeOptions(options,defaults);
assert(isscalar(options.hbar) && isfinite(options.hbar) && ...
    options.hbar > 0,'options.hbar must be positive.');
assert(isscalar(options.maximumExtendedDimension) && ...
    options.maximumExtendedDimension >= 1 && ...
    options.maximumExtendedDimension == ...
    fix(options.maximumExtendedDimension), ...
    'maximumExtendedDimension must be a positive integer.');

validationOptions = struct( ...
    'maximumOpticalDutyCycle',options.maximumOpticalDutyCycle, ...
    'requirePulsedOpticalSources',options.requirePulsedOpticalSources);
validation = driven_rttg.validateDriveSources( ...
    opticalSources,acousticSources,baseDimension,validationOptions);

relationOptions = struct( ...
    'relationSearchOrder',options.relationSearchOrder, ...
    'relationTolerance',options.relationTolerance, ...
    'maximumRelationCandidates',options.maximumRelationCandidates);
incommensurability = driven_rttg.checkDriveIncommensurability( ...
    opticalSources,acousticSources,relationOptions);
if options.requireNoDetectedRelations
    assert(~incommensurability.detectedRelation, ...
        ['A low-order combined space-time relation was detected: [%s]. ' ...
         'Change source vectors/frequencies or relax the diagnostic.'], ...
        sprintf(' %d',incommensurability.closestIntegerRelation));
end
if options.requireNoDetectedFrequencyRelations
    assert(~incommensurability.detectedFrequencyRelation, ...
        ['A low-order carrier-frequency relation was detected: [%s]. ' ...
         'The temporal phases are not independent within the search.'], ...
        sprintf(' %d',incommensurability.closestFrequencyRelation));
end

geometry = driven_rttg.getDriveSidebandGeometry( ...
    q,opticalSources,acousticSources,options.hbar);
if isnumeric(Hamiltonian)
    assert(all(vecnorm(geometry.waveVector,2,1) == 0), ...
        ['A fixed Hamiltonian matrix cannot represent nonzero drive ' ...
         'momentum transfer. Supply a function handle H(q).']);
end
extendedDimension = baseDimension*geometry.numberSidebands;
assert(extendedDimension <= options.maximumExtendedDimension, ...
    ['The extended dimension %d exceeds maximumExtendedDimension=%d. ' ...
     'Reduce source cutoffs or raise the explicit limit.'], ...
    extendedDimension,options.maximumExtendedDimension);

identity = speye(baseDimension);
diagonalBlocks = cell(geometry.numberSidebands,1);
for sideband = 1:geometry.numberSidebands
    Hq = evaluateHamiltonian(Hamiltonian,geometry.shiftedQ(:,sideband));
    assert(isequal(size(Hq),[baseDimension,baseDimension]), ...
        'The static Hamiltonian returned an unexpected dimension.');
    errorHermitian = norm(Hq-Hq','fro')/max(norm(Hq,'fro'),eps);
    assert(errorHermitian <= options.hermiticityTolerance, ...
        'The static Hamiltonian is not Hermitian at sideband %d.',sideband);
    diagonalBlocks{sideband} = sparse(Hq ...
        +geometry.sidebandEnergy(sideband)*identity);
end
staticMatrix = sparse(blkdiag(diagonalBlocks{:}));

sources = driven_rttg.combineDriveSources( ...
    opticalSources,acousticSources);
sourceMatrices = cell(numel(sources),1);
sidebandMap = makeSidebandMap(geometry.sidebandIndex);
for source = 1:numel(sources)
    if isempty(sources(source).couplingOperator)
        operator = identity;
    else
        operator = sparse(sources(source).couplingOperator);
    end
    positiveCoupling = (sources(source).amplitude/2) ...
        *exp(1i*sources(source).phase)*operator;
    adjacency = sparse(geometry.numberSidebands,geometry.numberSidebands);
    for sideband = 1:geometry.numberSidebands
        index = geometry.sidebandIndex(sideband,:);
        if index(source) >= sources(source).cutoff
            continue
        end
        neighbor = index;
        neighbor(source) = neighbor(source)+1;
        neighborSideband = sidebandMap(indexKey(neighbor));
        adjacency(sideband,neighborSideband) = 1;
    end
    forward = kron(adjacency,positiveCoupling);
    sourceMatrices{source} = sparse(forward+forward');
end

cache.staticMatrix = staticMatrix;
cache.sourceMatrices = sourceMatrices;
cache.geometry = geometry;
cache.opticalSources = opticalSources;
cache.acousticSources = acousticSources;
cache.sources = sources;
cache.baseDimension = baseDimension;
cache.extendedDimension = extendedDimension;
cache.centralSpatialRange = (geometry.centralSideband-1)*baseDimension ...
    +(1:baseDimension);
cache.validation = validation;
cache.incommensurability = incommensurability;
cache.options = options;

if options.verbose
    fprintf(['Prepared driven space-time Hamiltonian: %d spatial DoFs, ' ...
        '%d sidebands, extended dimension %d.\n'], ...
        baseDimension,geometry.numberSidebands,extendedDimension);
    fprintf('Closest searched source relation residual: %.3e.\n', ...
        incommensurability.minimumRelativeResidual);
    fprintf('%s\n',validation.warning);
end
end

function Hq = evaluateHamiltonian(Hamiltonian,q)
if isnumeric(Hamiltonian)
    Hq = Hamiltonian;
elseif isa(Hamiltonian,'function_handle')
    Hq = Hamiltonian(q);
else
    error('Hamiltonian must be a matrix or a function H(q).');
end
end

function map = makeSidebandMap(indices)
map = containers.Map('KeyType','char','ValueType','double');
for index = 1:size(indices,1)
    map(indexKey(indices(index,:))) = index;
end
end

function key = indexKey(index)
key = sprintf('%d,',index);
end
