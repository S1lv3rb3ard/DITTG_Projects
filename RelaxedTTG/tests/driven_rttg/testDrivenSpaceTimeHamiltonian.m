function tests = testDrivenSpaceTimeHamiltonian
%TESTDRIVENSPACETIMEHAMILTONIAN Locked indexing and pulse regression tests.
tests = functiontests(localfunctions);
end

function lockedDiagonalAndCouplingsAreCorrect(testCase)
baseDimension = 2;
Hfun = @(q) [q(1),1-0.2i;1+0.2i,q(2)];
q = [0.4;-0.3];
optical = driven_rttg.makeOpticalDriveSource( ...
    0.08,[0.11;-0.03],sqrt(2),0.2,1, ...
    struct('type','gaussianTrain','period',20,'fwhm',2,'center',1));
acoustic = driven_rttg.makeAcousticDriveSource( ...
    0.04,[-0.02;0.05],sqrt(3),-0.4,1);
options = struct( ...
    'hbar',0.7, ...
    'maximumOpticalDutyCycle',0.10, ...
    'relationSearchOrder',2, ...
    'maximumExtendedDimension',100, ...
    'verbose',false);
cache = driven_rttg.prepareDrivenSpaceTimeHamiltonian( ...
    Hfun,q,baseDimension,optical,acoustic,options);
[K,~,info] = driven_rttg.evaluateDrivenSpaceTimeHamiltonian(cache,1);

targetIndex = [1,-1];
sideband = find(all(cache.geometry.sidebandIndex == targetIndex,2));
range = (sideband-1)*baseDimension+(1:baseDimension);
shiftedQ = q-targetIndex(1)*optical.waveVector ...
    -targetIndex(2)*acoustic.waveVector;
energy = options.hbar*(targetIndex(1)*optical.frequency ...
    +targetIndex(2)*acoustic.frequency);
expectedDiagonal = Hfun(shiftedQ)+energy*eye(baseDimension);
verifyEqual(testCase,full(K(range,range)),expectedDiagonal,'AbsTol',1e-13);

zeroSideband = cache.geometry.centralSideband;
plusSideband = find(all(cache.geometry.sidebandIndex == [1,0],2));
zeroRange = (zeroSideband-1)*baseDimension+(1:baseDimension);
plusRange = (plusSideband-1)*baseDimension+(1:baseDimension);
expectedOptical = optical.amplitude/2*exp(1i*optical.phase) ...
    *eye(baseDimension);
verifyEqual(testCase,full(K(zeroRange,plusRange)), ...
    expectedOptical,'AbsTol',1e-13);
verifyLessThan(testCase,info.relativeHermiticityError,1e-14);

selector = driven_rttg.getDrivenCentralSelector( ...
    cache,speye(baseDimension));
verifyEqual(testCase,full(selector(zeroRange,:)),eye(baseDimension));
verifyEqual(testCase,nnz(selector),baseDimension);
end

function opticalPulseIsSuppressedBetweenShots(testCase)
optical = driven_rttg.makeOpticalDriveSource( ...
    1,[0.1;0.2],sqrt(2),0,1, ...
    struct('type','gaussianTrain','period',20,'fwhm',2,'center',0));
atPulse = driven_rttg.evaluateDriveEnvelope(optical,0);
betweenPulses = driven_rttg.evaluateDriveEnvelope(optical,10);
metrics = getOpticalPulseMetrics(optical);
verifyEqual(testCase,atPulse,1,'AbsTol',1e-14);
verifyLessThan(testCase,betweenPulses,1e-20);
verifyLessThan(testCase,metrics.equivalentIntensityDutyCycle,0.1);
end

function lowOrderRelationIsDetected(testCase)
optical = driven_rttg.makeOpticalDriveSource( ...
    1,[1;0],1,0,1, ...
    struct('type','gaussianTrain','period',20,'fwhm',2,'center',0));
acoustic = driven_rttg.makeAcousticDriveSource(1,[2;0],2,0,1);
report = driven_rttg.checkDriveIncommensurability( ...
    optical,acoustic,struct('relationSearchOrder',2));
verifyTrue(testCase,report.detectedRelation);
verifyEqual(testCase,report.spatialResidual,[0;0],'AbsTol',1e-14);
verifyEqual(testCase,report.frequencyResidual,0,'AbsTol',1e-14);
end

function multipleSourcesProduceExpectedSidebandCount(testCase)
envelope = struct( ...
    'type','gaussianTrain','period',30,'fwhm',2,'center',0);
optical(1) = driven_rttg.makeOpticalDriveSource( ...
    1,[1;0],sqrt(2),0,1,envelope);
optical(2) = driven_rttg.makeOpticalDriveSource( ...
    1,[0;1],sqrt(3),0,2,envelope);
acoustic(1) = driven_rttg.makeAcousticDriveSource( ...
    1,[1;1],sqrt(5),0,1);
geometry = driven_rttg.getDriveSidebandGeometry( ...
    [0;0],optical,acoustic,1);
verifyEqual(testCase,geometry.numberSidebands,3*5*3);
verifySize(testCase,geometry.sidebandIndex,[45,3]);
verifyEqual(testCase,numel(geometry.centralSideband),1);
end
