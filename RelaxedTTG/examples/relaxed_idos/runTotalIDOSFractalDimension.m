clearvars
close all

projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(projectRoot));

% Small demonstration of the total Green/Herglotz IDoS workflow.  The
% reciprocal-cell trapezoid is deliberately coarse and must be converged
% together with W, L, and the block-Lanczos depth in production.
a = 1.42*sqrt(3);
stack = getStack(a,[-0.3,0,0.3]);
W = 0.04;
L = 3;
DoF = getDoF(stack,W,L,'clean');
[qPoints,qWeights] = periodicCellTrapezoid(stack.B{2},3);

tA = [0,0.3302,0.23206,0.04969,-0.02499, ...
      0.00285,0.00204,-0.00014,-0.00029].';
tB = [-2.99251,-0.28983,0.02791,-0.00877,-0.01870, ...
      0.00621,-0.00256,-0.00018,-0.00033,-0.00264].';
shells = getInterpolatedIntralayerHoppingValues(stack,tA,tB,[]);
relaxationFields = makeToyTrilayerRelaxation(stack,0.01*a);

E = linspace(-4,4,4001);
eta = [0.08,0.04,0.02,0.01];

hamiltonianOptions.verbose = false;
hamiltonianOptions.intralayer = struct( ...
    'configurationGridSize',3,'useGPU',false,'verbose',false);
hamiltonianOptions.interlayer = struct( ...
    'realSpaceRadialOrder',12, ...
    'realSpaceAngularOrder',18, ...
    'realSpaceCutoff',6*a, ...
    'configurationGridSize',3, ...
    'minimumConfigurationGridSize',3, ...
    'spectatorModeCutoff',1, ...
    'modeBatchSize',4, ...
    'pairBatchSize',16, ...
    'useGPU',false, ...
    'verbose',false);

greenOptions.lanczos = struct( ...
    'maximumBlockSteps',40, ...
    'reorthogonalizationPasses',2, ...
    'breakdownTolerance',1e-12, ...
    'verbose',true);
greenOptions.energyChunkSize = 512;
greenOptions.storeLanczos = false;

options.hamiltonian = hamiltonianOptions;
options.green = greenOptions;
options.matrixCacheSize = 1;
options.verbose = true;

[greenResult,bridge] = relaxed_idos.computeRelaxedTotalIDOS( ...
    stack,DoF,shells,relaxationFields,qPoints,qWeights,E,eta,options); %#ok<ASGLU>

% Estimate a local measure dimension at E0 from IDoS increments.  The
% increment radius is kept four times larger than the Green broadening.
% A physical dimension requires repeating this calculation while increasing
% W, L, q resolution, and Lanczos depth as eta and the increment shrink.
scales = 4*eta(:);
dimensionOptions.evaluationEnergies = 0;
dimensionOptions.qOrders = [0,1,2];
dimensionOptions.window = [-1,1];
dimensionOptions.smoothingWidths = eta(:);
dimensionOptions.minimumScaleToSmoothingRatio = 4;
dimension = relaxed_idos.estimateFractalDimensionsFromIDOS( ...
    E,greenResult.totalIDOS,scales,dimensionOptions);

figure('Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile
plot(E,greenResult.totalIDOS,'LineWidth',1.2);
xlabel('energy (eV)');
ylabel('N_\eta(E)');
legend(compose('\\eta = %.3g eV',eta),'Location','best');
title('Green-operator integrated spectral measure');
grid on

nexttile
loglog(scales,dimension.local.mass(1,:),'o-','LineWidth',1.2);
xlabel('increment radius r (eV)');
ylabel('N_\eta(E_0+r)-N_\eta(E_0-r)');
title(sprintf('finite-scale local slope %.3f, R^2 %.3f', ...
    dimension.local.dimension(1),dimension.local.fitR2(1)));
grid on

fprintf('Estimated local dimension at E0=0: %.6f (R^2 %.4f)\n', ...
    dimension.local.dimension(1),dimension.local.fitR2(1));
fprintf('Generalized dimensions [D0,D1,D2]: [%s]\n', ...
    sprintf(' %.6f',dimension.generalized.dimension));
fprintf('%s\n',dimension.warning);
