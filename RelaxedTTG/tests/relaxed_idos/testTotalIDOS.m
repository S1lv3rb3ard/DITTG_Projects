function tests = testTotalIDOS
%TESTTOTALIDOS Exact small-matrix tests for the total Green IDoS.
tests = functiontests(localfunctions);
end

function blockGreenMatchesDirectResolvent(testCase)
rng(7);
dimension = 12;
randomMatrix = randn(dimension)+1i*randn(dimension);
H = 0.5*(randomMatrix+randomMatrix');
V = sparse([1,7],[1,2],1,dimension,2);
options = struct( ...
    'maximumBlockSteps',6, ...
    'reorthogonalizationPasses',2, ...
    'breakdownTolerance',1e-13, ...
    'verbose',false);
lanczos = rttg_common.blockLanczosHermitian(H,V,options);
z = [-0.4+0.2i,0.1+0.05i,0.8+0.3i];
[m,M] = rttg_common.evaluateBlockLanczosGreen(lanczos,z);

expectedM = complex(zeros(2,2,numel(z)));
expectedMTrace = complex(zeros(size(z)));
for index = 1:numel(z)
    expectedM(:,:,index) = V'*((H-z(index)*eye(dimension))\V);
    expectedMTrace(index) = trace(expectedM(:,:,index))/2;
end

verifyEqual(testCase,M,expectedM,'RelTol',2e-11,'AbsTol',2e-12);
verifyEqual(testCase,m,expectedMTrace,'RelTol',2e-11,'AbsTol',2e-12);
verifyGreaterThan(testCase,imag(m),zeros(size(m)));
verifyLessThan(testCase,lanczos.orthogonalityError,1e-11);
verifyGreaterThanOrEqual(testCase,lanczos.scalarWeights, ...
    zeros(size(lanczos.scalarWeights)));
verifyEqual(testCase,sum(lanczos.scalarWeights),1,'AbsTol',1e-13);
end

function greenIDOSMatchesExactPoissonFormula(testCase)
H = diag([-1.0,-0.25,0.5,1.25]);
V = speye(4);
E = linspace(-2,2,41);
eta = [0.1,0.03];
options.lanczos = struct( ...
    'maximumBlockSteps',2, ...
    'reorthogonalizationPasses',2, ...
    'verbose',false);
options.verbose = false;
result = relaxed_idos.computeTotalIDOS(H,[],1,V,E,eta,options);

lambda = diag(H);
expectedIDOS = zeros(numel(eta),numel(E));
expectedGreen = complex(zeros(size(expectedIDOS)));
for etaIndex = 1:numel(eta)
    expectedIDOS(etaIndex,:) = mean( ...
        0.5+atan((E-lambda)/eta(etaIndex))/pi,1);
    expectedGreen(etaIndex,:) = mean( ...
        1./(lambda-(E+1i*eta(etaIndex))),1);
end

verifyEqual(testCase,result.totalIDOS,expectedIDOS,'AbsTol',2e-13);
verifyEqual(testCase,result.green,expectedGreen,'AbsTol',2e-13);
verifyGreaterThanOrEqual(testCase,diff(result.totalIDOS,1,2), ...
    -1e-14*ones(numel(eta),numel(E)-1));
verifyEqual(testCase,sum(result.padeWeights),1,'AbsTol',1e-13);
end

function positiveQQuadratureProducesOnlyTotalIDOS(testCase)
H = {diag([-1,0]),diag([1,2])};
qPoints = [0,1;0,0];
qWeights = [1;3];
V = speye(2);
E = [-0.5,0.5,1.5];
eta = 0.1;
options.lanczos = struct('maximumBlockSteps',2,'verbose',false);
options.verbose = false;

result = relaxed_idos.computeTotalIDOS( ...
    H,qPoints,qWeights,V,E,eta,options);
localOne = mean(0.5+atan((E-diag(H{1}))/eta)/pi,1);
localTwo = mean(0.5+atan((E-diag(H{2}))/eta)/pi,1);
expected = 0.25*localOne+0.75*localTwo;

verifyEqual(testCase,result.totalIDOS,expected,'AbsTol',2e-13);
verifyFalse(testCase,isfield(result,'momentumLocalIDOS'));
verifyEqual(testCase,result.qWeights,[0.25;0.75],'AbsTol',eps);
end
