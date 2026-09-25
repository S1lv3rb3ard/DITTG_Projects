clearvars
close all

projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(projectRoot));

% A small, executable [-0.3,0,0.3]-degree example.  Increase W, L, the
% quadrature orders, P, and the number of points per segment for production.
a = 1.42*sqrt(3);
stack = getStack(a,[-1.4,0,2.8]);
W = 0.35;
L = 10;
DoF = getDoF(stack,W,L,'clean');

tA = [0,0.3302,0.23206,0.04969,-0.02499, ...
      0.00285,0.00204,-0.00014,-0.00029].';
tB = [-2.99251,-0.28983,0.02791,-0.00877,-0.01870, ...
      0.00621,-0.00256,-0.00018,-0.00033,-0.00264].';
shells = getInterpolatedIntralayerHoppingValues(stack,tA,tB,[]);

% Replace this software-test field with the mechanically minimized fields
% when comparing with physical spectra.
relaxationFields = makeToyTrilayerRelaxation(stack,0.01*a);

pointsPerSegment = 8;
linecut = rttg_common.getLinecut(stack.K,pointsPerSegment);
E = linspace(-0.2, 0.2, 5e3);

hamiltonianOptions.verbose = false;
hamiltonianOptions.intralayer = struct( ...
    'configurationGridSize', 10, ...
    'couplingInnerRadius', 6.00, ...
    'couplingOuterRadius', 6.25, ...
    'useGPU', true, ...
    'verbose', true);
hamiltonianOptions.interlayer = struct( ...
    'realSpaceRadialOrder', 80, ...
    'realSpaceAngularOrder', 124, ...
    'realSpaceCutoff', 8*a, ...
    'configurationGridSize', 15, ...
    'minimumConfigurationGridSize', 3, ...
    'spectatorModeCutoff', inf, ...
    'modeBatchSize', 4, ...
    'pairBatchSize', 32, ...
    'useGPU', true, ...
    'verbose', true);

% P is the maximum Chebyshev order.  P=4000 is large enough to make this a
% meaningful spectral diagnostic; use a smaller value only for code-path
% smoke tests.
ldosOptions.P = 5e3;
ldosOptions.mode = 'gpu';
ldosOptions.scale = []; % rigorous automatic bound over this line cut
ldosOptions.scaleSafetyFactor = 1.25;
% This example is small enough to retain every assembled line-cut matrix,
% avoiding a second H(q) assembly during the KPM pass.
ldosOptions.matrixCacheSize = size(linecut.q,2);
ldosOptions.W = W;
ldosOptions.L = L;
ldosOptions.hamiltonian = hamiltonianOptions;
ldosOptions.verbose = true;

[result,Job,bridge] = relaxed_dos.computeMomentumLocalDOSLinecut( ...
    stack,DoF,shells,relaxationFields,linecut,E,ldosOptions);

fprintf('Estimated central KPM energy resolution: %.4g eV\n', ...
    result.estimatedCentralEnergyResolution);

plotOptions.title = ...
    'Relaxed momentum LDoS, angles [-0.3, 0, 0.3] degrees';
relaxed_dos.plotMomentumLocalDOSLinecut( ...
    linecut,E,result.LDoS,plotOptions);
