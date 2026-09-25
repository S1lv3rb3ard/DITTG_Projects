function tests = testMomentumLocalIDOS
%TESTMOMENTUMLOCALIDOS Preserve q dependence in the Green IDoS workflow.
tests = functiontests(localfunctions);
end

function linecutMatchesExactPoissonFormula(testCase)
H = {diag([-1,0.25,1.5]),diag([-0.5,0.75,2.0])};
qPoints = [0,1;0,0];
V = speye(3);
E = linspace(-2,2.5,31);
eta = [0.2,0.05];
options.lanczos = struct( ...
    'maximumBlockSteps',2, ...
    'reorthogonalizationPasses',2, ...
    'verbose',false);
options.verbose = false;

result = relaxed_idos.computeMomentumLocalIDOS( ...
    H,qPoints,V,E,eta,options);
expected = zeros(2,numel(E),numel(eta));
for qIndex = 1:2
    lambda = diag(H{qIndex});
    for etaIndex = 1:numel(eta)
        expected(qIndex,:,etaIndex) = mean( ...
            0.5+atan((E-lambda)/eta(etaIndex))/pi,1);
    end
end

verifyEqual(testCase,result.momentumLocalIDOS,expected,'AbsTol',2e-13);
verifySize(testCase,result.momentumLocalIDOS,[2,numel(E),numel(eta)]);
verifyGreaterThanOrEqual(testCase,diff(result.momentumLocalIDOS,1,2), ...
    -1e-14*ones(2,numel(E)-1,numel(eta)));
verifyNotEqual(testCase,result.momentumLocalIDOS(1,:,:), ...
    result.momentumLocalIDOS(2,:,:));
verifyFalse(testCase,isfield(result,'totalIDOS'));
end

function surfacePlotUsesSelectedEtaSlice(testCase)
linecut.q = [0,1;0,0];
linecut.kPath = [0,1];
linecut.tickLocs = [0,1];
linecut.labels = {'$K_1$','$K_2$'};
E = [-1,0,1];
eta = [0.2,0.05];
IDOS = zeros(2,3,2);
IDOS(:,:,1) = [0.1,0.4,0.9;0.2,0.5,0.8];
IDOS(:,:,2) = [0.05,0.35,0.95;0.15,0.55,0.85];

options = struct('figureVisible','off','etaIndex',2);
[fig,~,surfaceHandle] = relaxed_idos.plotMomentumLocalIDOSLinecut( ...
    linecut,E,IDOS,eta,options);
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
verifyEqual(testCase,surfaceHandle.ZData,IDOS(:,:,2).');
verifyEqual(testCase,surfaceHandle.CData,IDOS(:,:,2).');
end
