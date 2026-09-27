function tests = testRelaxationConventions
%TESTRELAXATIONCONVENTIONS Julia loader, disregistry points, intralayer gauge.
%
% The fixture tests/data/triG_*_1.50_2.00_6.jld is a converged N = 6
% relaxation for theta12 = 1.5, theta23 = 2.0 degrees in the exact binary
% layout written by example2.jl. With the Julia rotation convention this
% corresponds to getStack(a,[-1.5,0,-2.0]).
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.dataFolder = fullfile(root,'data');
testCase.TestData.a = 1.42*sqrt(3);
testCase.TestData.stack = getStack(testCase.TestData.a,[-1.5,0,-2.0]);
testCase.TestData.minimizer = fullfile(testCase.TestData.dataFolder, ...
    'triG_minimizer_1.50_2.00_6.jld');
testCase.TestData.data = fullfile(testCase.TestData.dataFolder, ...
    'triG_data_1.50_2.00_6.jld');
end

function testLoaderReproducesHullValues(testCase)
% At the configuration of a Julia hull node gamma (fractional coordinates of
% the shift b = R - x), RelaxedTTG coordinates are x = -Phi*tE_t*gamma and
% the field must equal Phi*u_Julia(node) exactly.
stack = testCase.TestData.stack;
[fields,info] = loadJuliaTrilayerRelaxation(testCase.TestData.minimizer, ...
    testCase.TestData.data,stack,struct('relativeTolerance',0,'verbose',false));
verifyEqual(testCase,info.frameAngleDegrees,30,'AbsTol',1e-10);
julia = readJuliaRelaxationFiles(testCase.TestData.minimizer, ...
    testCase.TestData.data);
N = julia.N;
jrot = @(t) [cos(t),sin(t);-sin(t),cos(t)];
thetaJulia = [julia.theta1,0,julia.theta3];
juliaPairs = {[2,3],[1,3],[2,1]};
rng(3);
for layer = 1:3
    others = setdiff(1:3,layer,'stable');
    nodes = randi(N,4,25)-1;
    x = cell(1,3);
    for slot = 1:2
        relative = others(slot);
        pairIndex = find(juliaPairs{layer} == relative);
        gamma = nodes(2*pairIndex-1:2*pairIndex,:)/N;
        x{relative} = -info.Phi*(jrot(thetaJulia(relative))*julia.E)*gamma;
    end
    computed = fields{layer}(x{others(1)},x{others(2)});
    expected = zeros(2,size(nodes,2));
    for p = 1:size(nodes,2)
        index = num2cell(nodes(:,p)+1);
        expected(:,p) = info.Phi*julia.u(:,index{:},layer);
    end
    verifyEqual(testCase,computed,expected,'AbsTol',1e-12);
end
end

function testLoaderRejectsMismatchedStack(testCase)
stack = getStack(testCase.TestData.a,[1.5,0,2.0]);   % wrong rotation sense
verifyError(testCase,@() loadJuliaTrilayerRelaxation( ...
    testCase.TestData.minimizer,testCase.TestData.data,stack, ...
    struct('verbose',false)),'loadJuliaTrilayerRelaxation:geometryMismatch');
end

function testGetRelaxationFieldsFindsJuliaOutput(testCase)
[~,info] = getRelaxationFields(testCase.TestData.stack, ...
    testCase.TestData.dataFolder,struct('loader',struct('verbose',false)));
verifyEqual(testCase,info.source,'julia');
verifyEqual(testCase,info.N,6);
end

function testSublatticeBondsCarryOnlyPhysicalStrain(testCase)
% With disregistry-consistent evaluation, the shortest bond of each channel
% changes by the continuum strain times its length. For this 1.5/2.0 degree
% fixture, max|u| ~ 0.08 A over a moire length ~ 94 A gives a strain of
% ~5e-3 and changes of ~0.01 A on the 2.46 A lattice bonds (AA, BB) and
% less on the 1.42 A A-B bonds. Evaluating at T_j(R+tau) would change the
% bonds by O(max|u|) ~ 0.08 A, so 0.03 A separates the two conventions.
stack = testCase.TestData.stack;
fields = loadJuliaTrilayerRelaxation(testCase.TestData.minimizer, ...
    testCase.TestData.data,stack,struct('verbose',false));
tA = [0,0.3302,0.23206].';
tB = [-2.99251,-0.28983,0.02791].';
shells = getInterpolatedIntralayerHoppingValues(stack,tA,tB,[]);
[s1,s2,s3,s4] = ndgrid((0:4)/5);
for layer = 1:3
    others = setdiff(1:3,layer,'stable');
    bk = stack.A{others(1)}*[s1(:),s2(:)].';
    bl = stack.A{others(2)}*[s3(:),s4(:)].';
    channels = sampleRelaxedIntralayerChannels(stack,shells, ...
        fields{layer},layer,bk,bl,struct());
    for c = 1:4
        length0 = vecnorm(channels(c).bond0,2,1);
        nearest = abs(length0-min(length0(length0 > 1e-9))) < 1e-9;
        change = abs(channels(c).relaxedLength(:,nearest)- ...
            vecnorm(channels(c).bond0(:,nearest),2,1));
        verifyLessThan(testCase,max(change,[],'all'),3e-2, ...
            sprintf('layer %d channel %s',layer,channels(c).name));
    end
end
end

function testUnrelaxedIntralayerUsesOneBlochGauge(testCase)
% With zero relaxation the central 2x2 intralayer block of layer j must be
%   H_AA(q) = sum_d t(|d|) exp(-i q.d),
%   H_AB(q) = sum_{r in d+tau} t(|r|) exp(+i q.r)  (row A, column B),
% i.e. the same Bloch convention (orbital offsets included) as the
% interlayer blocks. A mixed gauge makes the Dirac cone anisotropic.
a = testCase.TestData.a;
stack = getStack(a,[-1.4,0,2.8]);
tA = [0,0.3302,0.23206,0.04969].';
tB = [-2.99251,-0.28983,0.02791,-0.00877].';
shells = getInterpolatedIntralayerHoppingValues(stack,tA,tB,[]);
DoF = getDoF(stack,0.25,3,'clean');
zeroField = @(xk,xl) zeros(size(xk),'like',xk);
cache = prepareRelaxedIntralayerHamiltonian(stack,DoF,shells, ...
    {zeroField,zeroField,zeroField}, ...
    struct('configurationGridSize',2,'verbose',false));
rng(5);
for layer = 1:3
    central = find(DoF(:,9) == layer & all(DoF(:,3:8) == 0,2));
    assumeNotEmpty(testCase,central);
    shellMap = stack.A{layer}/stack.A{2};
    RA = shellMap*shells.matA;
    RB = shellMap*shells.matB;
    for trial = 1:3
        q = stack.K(:,layer)+0.05*randn(2,1);
        H = evaluateRelaxedIntralayerHamiltonian(cache,q);
        rows = 2*central-1+(0:1);
        block = full(H(rows,rows));
        expectedAA = sum(shells.intraAA(vecnorm(RA,2,1)).*exp(-1i*q.'*RA));
        expectedAB = sum(shells.intraAB(vecnorm(RB,2,1)).*exp(1i*q.'*RB));
        verifyEqual(testCase,block(1,1),expectedAA,'AbsTol',1e-12);
        verifyEqual(testCase,block(2,2),expectedAA,'AbsTol',1e-12);
        verifyEqual(testCase,block(1,2),expectedAB,'AbsTol',1e-12);
        verifyEqual(testCase,block(2,1),conj(expectedAB),'AbsTol',1e-12);
    end
end
end
