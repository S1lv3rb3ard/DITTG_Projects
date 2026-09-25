function tests = testMomentumLocalDOSLinecut
%TESTMOMENTUMLOCALDOSLINECUT End-to-end line-cut and plotting tests.
tests = functiontests(localfunctions);
end

function momentumLDoSSelectorHasSixOrbitals(testCase)
stack = getStack(1.42*sqrt(3),[-1.4,0,2.8]);
DoF = getDoF(stack,0.25,3,'clean');
X = getCenterBasis(DoF,1:6);
verifySize(testCase,X,[2*size(DoF,1),6]);
verifyEqual(testCase,full(sum(X,1)),ones(1,6));
verifyEqual(testCase,full(sum(X,2)),double(sum(X,2) > 0));
end

function highSymmetryPathAndPlotAreConsistent(testCase)
stack = getStack(1.42*sqrt(3),[-0.3,0,0.3]);
pointsPerSegment = 2;
linecut = rttg_common.getLinecut(stack.K,pointsPerSegment);

verifySize(testCase,linecut.q,[2,5*pointsPerSegment+1]);
verifyEqual(testCase,linecut.q(:,1),stack.K(:,1),'AbsTol',1e-14);
verifyEqual(testCase,linecut.q(:,end),stack.K(:,3),'AbsTol',1e-14);
verifyGreaterThan(testCase,diff(linecut.kPath),zeros(1,5*pointsPerSegment));
verifyEqual(testCase,linecut.kPath([1,end]), ...
    linecut.tickLocs([1,end]),'AbsTol',1e-14);

E = linspace(-0.2,0.2,9);
LDoS = reshape(1:(numel(linecut.kPath)*numel(E)), ...
    numel(linecut.kPath),numel(E));
[fig,ax,imageHandle] = relaxed_dos.plotMomentumLocalDOSLinecut( ...
    linecut,E,LDoS,struct('figureVisible','off'));
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
verifyEqual(testCase,imageHandle.CData,LDoS.');
verifyEqual(testCase,ax.XTick,linecut.tickLocs,'AbsTol',1e-14);
end

function relaxedLinecutProducesFiniteMomentumLDoS(testCase)
a = 1.42*sqrt(3);
stack = getStack(a,[-1.4,0,2.8]);
W = 0.10;
L = 3;
DoF = getDoF(stack,W,L,'clean');
tA = [0,0.3302,0.23206,0.04969,-0.02499, ...
      0.00285,0.00204,-0.00014,-0.00029].';
tB = [-2.99251,-0.28983,0.02791,-0.00877,-0.01870, ...
      0.00621,-0.00256,-0.00018,-0.00033,-0.00264].';
shells = getInterpolatedIntralayerHoppingValues(stack,tA,tB,[]);
fields = makeToyTrilayerRelaxation(stack,0.002*a);
linecut = rttg_common.getLinecut(stack.K,1);
E = linspace(-0.2,0.2,11);

hamiltonianOptions.verbose = false;
hamiltonianOptions.intralayer = struct( ...
    'configurationGridSize',2,'useGPU',false,'verbose',false);
hamiltonianOptions.interlayer = struct( ...
    'realSpaceRadialOrder',4, ...
    'realSpaceAngularOrder',6, ...
    'realSpaceCutoff',4*a, ...
    'configurationGridSize',3, ...
    'minimumConfigurationGridSize',3, ...
    'spectatorModeCutoff',0, ...
    'modeBatchSize',2, ...
    'pairBatchSize',8, ...
    'useGPU',false, ...
    'verbose',false);
options = struct( ...
    'P',16, ...
    'mode','cpu', ...
    'matrixCacheSize',size(linecut.q,2), ...
    'W',W, ...
    'L',L, ...
    'hamiltonian',hamiltonianOptions, ...
    'verbose',false);

[result,Job,bridge] = ...
    relaxed_dos.computeMomentumLocalDOSLinecut( ...
    stack,DoF,shells,fields,linecut,E,options);

verifySize(testCase,result.LDoS,[size(linecut.q,2),numel(E)]);
verifyTrue(testCase,all(isfinite(result.LDoS),'all'));
verifyLessThan(testCase,result.scale*result.spectralRadiusBound,1);
verifyEqual(testCase,result.estimatedCentralEnergyResolution, ...
    pi/((options.P+1)*result.scale),'RelTol',1e-14);
verifySize(testCase,result.projector,[bridge.dimension,6]);
verifyEqual(testCase,Job.list.Q,linecut.q.','AbsTol',0);
verifyFalse(testCase,ismember('cost',Job.list.Properties.VariableNames));
end
