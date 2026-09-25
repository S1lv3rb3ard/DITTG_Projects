function [result,Job,bridge] = computeMomentumLocalDOSLinecut( ...
    stack,DoF,shellsReference,relaxationFields,linecut,E,options)
%COMPUTEMOMENTUMLOCALDOSLINECUT Compute the paper's LDoS along a path.
%
% [result,Job,bridge] = relaxed_dos.computeMomentumLocalDOSLinecut(...)
% constructs the cached relaxed Hamiltonian on linecut.q, selects the six
% (G=0,j,alpha) orbitals, computes their KPM moments, and reconstructs
%
%   (1/6) sum_{j,alpha} [delta_epsilon(E-H(q))]_(0,j,alpha),(0,j,alpha).
%
% This is the momentum LDoS of equations (3.19), (3.20), and (4.51), not a
% real-space LDOS.  The relaxed interlayer transform used by H(q) is the
% direct continuous quadrature implemented by prepareRelaxedHamiltonian;
% no FFT is used.
if nargin < 7
    options = struct;
end
defaults = struct( ...
    'P',512, ...
    'mode','cpu', ...
    'scale',[], ...
    'validateProvidedScale',false, ...
    'scaleSafetyFactor',1.25, ...
    'matrixCacheSize',1, ...
    'W',NaN, ...
    'L',NaN, ...
    'hamiltonian',struct, ...
    'verbose',true);
options = mergeOptions(options,defaults);

validateInputs(linecut,E,options);
E = E(:).';
numberQ = size(linecut.q,2);

JobList = table;
JobList.t = ones(numberQ,1);
JobList.X = ones(numberQ,1);
JobList.N = ones(numberQ,1);
JobList.Q = linecut.q.';
JobList.cost = zeros(numberQ,1);
JobList.model = repmat("relaxed",numberQ,1);

Job.list = JobList;
Job.theta = stack.theta;
Job.W = options.W;
Job.L = options.L;
Job.P = options.P;
Job.H = cell(1,1);
Job.X = cell(1,1);
Job.X{1} = {getCenterBasis(DoF,1:6)};

bridgeOptions = struct( ...
    'hamiltonian',options.hamiltonian, ...
    'matrixCacheSize',options.matrixCacheSize, ...
    'validateSelectors',true, ...
    'setMomentumLDoSSelector',true, ...
    'momentumLDoSSelectorIndex',1, ...
    'verbose',options.verbose);
[Job,bridge] = relaxed_dos.attachRelaxedHamiltonianToDoS( ...
    Job,1,stack,DoF,shellsReference,relaxationFields,bridgeOptions);

% The infinity norm is a rigorous spectral-radius bound for Hermitian H.
% Automatic scaling therefore scans the whole path.  A supplied scale can
% skip that extra pass; set validateProvidedScale=true to check it here.
scanAllMatrices = isempty(options.scale) || options.validateProvidedScale;
spectralBounds = NaN(numberQ,1);
cost = zeros(numberQ,1);
if scanAllMatrices
    for qIndex = 1:numberQ
        Hq = bridge.H(linecut.q(:,qIndex));
        spectralBounds(qIndex) = norm(Hq,inf);
        cost(qIndex) = estimateKPMCost(Hq,6,options.P);
    end
    maximumBound = max(spectralBounds);
else
    Hq = bridge.H(linecut.q(:,1));
    spectralBounds(1) = norm(Hq,inf);
    cost(:) = estimateKPMCost(Hq,6,options.P);
    maximumBound = NaN;
end

if isempty(options.scale)
    assert(maximumBound > 0, ...
        ['The Hamiltonian is zero on the complete line cut; automatic ' ...
         'scaling is undefined.']);
    scale = 1/(options.scaleSafetyFactor*maximumBound);
else
    scale = options.scale;
end
assert(isscalar(scale) && scale > 0 && isfinite(scale), ...
    'options.scale must be a positive scalar.');
if scanAllMatrices
    assert(scale*maximumBound < 1, ...
        ['The requested scale does not map every sampled H(q) into ' ...
         '(-1,1). Use scale < %.16g.'],1/maximumBound);
end
assert(max(abs(scale*E)) < 1, ...
    ['The requested energy grid lies outside the scaled KPM interval. ' ...
     'Require max(abs(E)) < %.16g.'],1/scale);

Job.scale = scale;
Job.list.cost = cost;
Job = relaxed_dos.computeMoments(Job,options.mode);
Job = relaxed_dos.computeMomentumLocalDOS(Job,E,options.mode);

result.energy = E;
result.q = linecut.q;
result.kPath = linecut.kPath;
result.tickLocs = linecut.tickLocs;
result.labels = linecut.labels;
result.LDoS = Job.list.DoS;
result.scale = scale;
result.estimatedCentralEnergyResolution = pi/((options.P+1)*scale);
result.spectralRadiusBound = maximumBound;
result.sampledSpectralRadiusBounds = spectralBounds;
result.projector = bridge.momentumLDoSBasis;

if options.verbose
    fprintf(['Computed relaxed momentum LDoS at %d q-points and %d ' ...
        'energies (P=%d, scale=%.6g, estimated central resolution ' ...
        '%.4g eV).\n'],numberQ,numel(E),options.P,scale, ...
        result.estimatedCentralEnergyResolution);
end
end

function bytes = estimateKPMCost(H,numberSelected,P)
% Complex sparse storage plus the block-Chebyshev work arrays and moments.
dimension = size(H,1);
sparseBytes = 24*nnz(H)+8*(dimension+1);
workBytes = 5*16*dimension*numberSelected;
momentBytes = 8*(P+1);
bytes = 1.25*(2*sparseBytes+workBytes+momentBytes);
end

function validateInputs(linecut,E,options)
assert(isstruct(linecut) && all(isfield(linecut, ...
    {'q','kPath','tickLocs','labels'})), ...
    'linecut must be produced by getLinecut.');
assert(size(linecut.q,1) == 2,'linecut.q must be 2-by-N.');
assert(numel(linecut.kPath) == size(linecut.q,2), ...
    'linecut.kPath and linecut.q must contain the same number of points.');
assert(isvector(E) && ~isempty(E) && all(isfinite(E)), ...
    'E must be a finite, nonempty vector.');
assert(options.P >= 1 && options.P == fix(options.P), ...
    'options.P must be a positive integer.');
assert(any(strcmpi(options.mode,{'cpu','gpu'})), ...
    'options.mode must be cpu or gpu.');
assert(isscalar(options.validateProvidedScale) && ...
    (islogical(options.validateProvidedScale) || ...
     (isnumeric(options.validateProvidedScale) && ...
      any(options.validateProvidedScale == [0,1]))), ...
    'options.validateProvidedScale must be a logical scalar.');
assert(isscalar(options.scaleSafetyFactor) && ...
    options.scaleSafetyFactor > 1, ...
    'options.scaleSafetyFactor must exceed one.');
assert(options.matrixCacheSize >= 0 && ...
    options.matrixCacheSize == fix(options.matrixCacheSize), ...
    'options.matrixCacheSize must be a nonnegative integer.');
end
