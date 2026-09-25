clearvars
close all

projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(projectRoot));

% Units in this example are eV, Angstrom, and fs. Frequencies are angular
% frequencies in rad/fs, so hbar converts sideband frequencies to eV.
hbar = 0.6582119569; % eV fs
a = 1.42*sqrt(3);
stack = getStack(a,[-0.3,0,0.3]);
W = 0.02;
L = 3;
DoF = getDoF(stack,W,L,'clean');
q = stack.K(:,2);

tA = [0,0.3302,0.23206,0.04969,-0.02499, ...
      0.00285,0.00204,-0.00014,-0.00029].';
tB = [-2.99251,-0.28983,0.02791,-0.00877,-0.01870, ...
      0.00621,-0.00256,-0.00018,-0.00033,-0.00264].';
shells = getInterpolatedIntralayerHoppingValues(stack,tA,tB,[]);
relaxationFields = makeToyTrilayerRelaxation(stack,0.005*a);

% Two pulsed optical carriers.  The periods, carrier frequencies, and wave
% vectors are intentionally different.  The amplitudes are Hamiltonian
% coupling parameters, not calibrated laser intensities.
optical(1) = driven_rttg.makeOpticalDriveSource( ...
    0.020,[0.013;0.004],0.40/hbar,0.17,1, ...
    struct('type','gaussianTrain','period',80,'fwhm',8,'center',20));
optical(2) = driven_rttg.makeOpticalDriveSource( ...
    0.014,[-0.007;0.016],sqrt(2)*0.27/hbar,-0.31,1, ...
    struct('type','sin2Train','period',90*sqrt(2), ...
           'duration',12,'center',35));

% Two coherent acoustic deformation-potential sources.
acoustic(1) = driven_rttg.makeAcousticDriveSource( ...
    0.004,[0.0021;0.0007],sqrt(3)*0.004/hbar,0.23,1);
acoustic(2) = driven_rttg.makeAcousticDriveSource( ...
    0.003,[-0.0012;0.0026],sqrt(5)*0.003/hbar,-0.11,1);

hamiltonianOptions.verbose = false;
hamiltonianOptions.intralayer = struct( ...
    'configurationGridSize',3,'useGPU',false,'verbose',false);
hamiltonianOptions.interlayer = struct( ...
    'realSpaceRadialOrder',10, ...
    'realSpaceAngularOrder',16, ...
    'realSpaceCutoff',5*a, ...
    'configurationGridSize',3, ...
    'minimumConfigurationGridSize',3, ...
    'spectatorModeCutoff',1, ...
    'modeBatchSize',4, ...
    'pairBatchSize',16, ...
    'useGPU',false, ...
    'verbose',false);

drivingOptions.hbar = hbar;
drivingOptions.maximumExtendedDimension = 10000;
drivingOptions.maximumOpticalDutyCycle = 0.10;
drivingOptions.requirePulsedOpticalSources = true;
drivingOptions.requireNoDetectedRelations = true;
drivingOptions.relationSearchOrder = 3;
drivingOptions.relationTolerance = 1e-10;
drivingOptions.verbose = true;

options.hamiltonian = hamiltonianOptions;
options.driving = drivingOptions;
options.matrixCacheSize = 1;

[driveCache,bridge] = ...
    driven_rttg.prepareRelaxedDrivenSpaceTimeHamiltonian( ...
    stack,DoF,shells,relaxationFields,q,optical,acoustic,options); %#ok<ASGLU>

% The carrier oscillations are represented by the sideband indices.  This
% time variable controls the slow optical pulse envelopes.
envelopeTime = 35;
[K,parts,info] = driven_rttg.evaluateDrivenSpaceTimeHamiltonian( ...
    driveCache,envelopeTime); %#ok<ASGLU>

fprintf('Driven dimension: %d\n',info.dimension);
fprintf('Driven nonzeros: %d\n',info.numNonzeros);
fprintf('Relative Hermiticity error: %.3e\n', ...
    info.relativeHermiticityError);
fprintf('Envelope values at t=%.3f fs: [%s]\n', ...
    envelopeTime,sprintf(' %.5f',info.envelopeValue));

time = linspace(0,240,1500);
figure('Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile
hold on
for source = 1:numel(optical)
    plot(time,driven_rttg.evaluateDriveEnvelope( ...
        optical(source),time), ...
        'LineWidth',1.2,'DisplayName',sprintf('optical %d',source));
end
xlabel('envelope time (fs)');
ylabel('field envelope');
ylim([-0.02,1.05]);
legend('Location','best');
grid on
title('Pulsed optical sources');

nexttile
spy(K);
title(sprintf('K_q(t), dimension %d',size(K,1)));

% Evaluate the note's central-sector observable.  These modest Lanczos and
% sideband cutoffs make this an executable demonstration, not a converged
% production spectrum.
observableEnergy = linspace(-0.8,0.8,601);
observableEta = 0.04;
observableOptions.lanczos = struct( ...
    'maximumBlockSteps',16, ...
    'reorthogonalizationPasses',2, ...
    'verbose',false);
observableOptions.computeOpticalAcousticContrast = true;
observableOptions.verbose = true;
observable = driven_rttg.computeDrivenSpectralObservables( ...
    driveCache,observableEnergy,observableEta,envelopeTime, ...
    observableOptions);

spectralFunction = squeeze( ...
    observable.combined.totalSpectralFunction(1,1,:));
integratedMeasure = squeeze( ...
    observable.combined.totalIntegratedMeasure(1,1,:));
contrast = squeeze(observable.contrast.totalSpectralFunction(1,1,:));

figure('Color','w');
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile
plot(observableEnergy,spectralFunction,'LineWidth',1.2);
xlabel('E (eV)');
ylabel('A_{00}(q,E)');
grid on
title('Central-sector spectrum');

nexttile
plot(observableEnergy,integratedMeasure,'LineWidth',1.2);
xlabel('E (eV)');
ylabel('N_{00}(q,E)');
ylim([0,1]);
grid on
title('Integrated measure');

nexttile
plot(observableEnergy,contrast,'LineWidth',1.2);
xlabel('E (eV)');
ylabel('C_{\gamma ac}(q,E)');
grid on
title('Optical-acoustic contrast');

boundaryDiagnostic = observable.diagnostics{1}.combined;
fprintf('Central-measure-weighted sideband-boundary weight: %.3e\n', ...
    boundaryDiagnostic.centralMeasureWeightedBoundaryWeight);

% The pulse-duty guard does not predict a damage threshold.  Experimental
% peak intensity and repetition rate must be checked using a device-specific
% absorption and heat-flow model before applying a laboratory drive.
