function tests = testInterlayerCompression
%TESTINTERLAYERCOMPRESSION Uniform compression of the interlayer hopping.
tests = functiontests(localfunctions);
end

function carrFitsAtZeroAndTenPercent(testCase)
p0 = getInterlayerHoppingParameters('carr2018',0);
verifyEqual(testCase,[p0.lambda0,p0.xi0,p0.kappa0,p0.lambda3,p0.xi3, ...
    p0.x3,p0.lambda6,p0.xi6,p0.x6,p0.kappa6], ...
    [0.310,1.750,1.990,-0.068,3.286,0.500,-0.008,2.272,1.217,1.562], ...
    'AbsTol',1e-15);
% Compression strengthens the coupling: lambda0(-0.1) = 0.310+0.1882+0.07741.
p = getInterlayerHoppingParameters('carr2018',-0.1);
verifyEqual(testCase,p.lambda0,0.310+0.1882+0.07741,'AbsTol',1e-12);
verifyGreaterThan(testCase,p.lambda0,p0.lambda0);
end

function pressureLawRoundTrip(testCase)
epsilon = compressionFromPressure(9.2);
verifyEqual(testCase,epsilon,-0.1,'AbsTol',2e-3);
verifyEqual(testCase,5.73*(exp(-9.54*epsilon)-1),9.2,'AbsTol',1e-12);
verifyEqual(testCase,compressionFromPressure(0),0);
end

function extrapolationWarns(testCase)
verifyWarning(testCase,@() getInterlayerHoppingParameters('carr2018',0.2), ...
    'getInterlayerHoppingParameters:extrapolation');
verifyWarning(testCase,@() getInterlayerHoppingParameters('carr2018',-0.25), ...
    'getInterlayerHoppingParameters:extrapolation');
verifyError(testCase,@() getInterlayerHoppingParameters('fangKaxiras2016',-0.1), ...
    'MATLAB:assertion:failed');
end

function legacyModelReproducesPreviousHopping(testCase)
stack = getStack(1.42*sqrt(3),[-1.4,0,2.8]);
stack = setInterlayerCompression(stack,0,'fangKaxiras2016');
rng(7);
r = 3*randn(2,200);
for alpha = 1:2
    for beta = 1:2
        h = realSpaceInterlayerHopping(r,stack,1,2,alpha,beta);
        verifyEqual(testCase,h,legacyHopping(r,stack,1,2,alpha,beta), ...
            'AbsTol',1e-15);
    end
end
end

function compressionIsPerInterface(testCase)
stack = getStack(1.42*sqrt(3),[-1.4,0,2.8]);
stack = setInterlayerCompression(stack,[-0.1,0]);
r = [0.4;0.1];
compressed = realSpaceInterlayerHopping(r,stack,1,2,1,1);
equilibrium = realSpaceInterlayerHopping(r,stack,2,3,1,1);
reference = realSpaceInterlayerHopping(r,getStack(1.42*sqrt(3), ...
    [-1.4,0,2.8]),2,3,1,1);
verifyGreaterThan(testCase,abs(compressed),abs(equilibrium));
verifyEqual(testCase,equilibrium,reference,'AbsTol',1e-15);
verifyError(testCase,@() realSpaceInterlayerHopping(r,stack,1,3,1,1), ...
    'realSpaceInterlayerHopping:unsupportedPair');
end

function expansionContinuationIsSmoothAndDecaying(testCase)
% C^1 at eps = 0, amplitudes decrease monotonically with distance over the
% whole supported range, and the array form matches the scalar form.
names = {'lambda0','xi0','kappa0','lambda3','xi3','x3', ...
    'lambda6','xi6','x6','kappa6'};
h = 1e-7;
below = getInterlayerHoppingParameters('carr2018',-h);
at = getInterlayerHoppingParameters('carr2018',0);
above = getInterlayerHoppingParameters('carr2018',h);
for n = 1:numel(names)
    verifyEqual(testCase,above.(names{n}),below.(names{n}),'AbsTol',1e-6);
    slopeBelow = (at.(names{n})-below.(names{n}))/h;
    slopeAbove = (above.(names{n})-at.(names{n}))/h;
    verifyEqual(testCase,slopeAbove,slopeBelow,'AbsTol',1e-4);
end
epsilon = linspace(-0.2,0.15,351);
p = getInterlayerHoppingParameters('carr2018',epsilon);
verifyTrue(testCase,all(diff(p.lambda0) < 0));
verifyTrue(testCase,all(diff(abs(p.lambda3)) < 0));
verifyTrue(testCase,all(diff(abs(p.lambda6)) < 0));
single = getInterlayerHoppingParameters('carr2018',epsilon(300));
verifyEqual(testCase,p.lambda0(300),single.lambda0,'AbsTol',1e-15);
end

function localDistanceOverridesInterface(testCase)
stack = getStack(1.42*sqrt(3),[-1.4,0,2.8]);
r = [0.3,1.2,-0.7;0.1,-0.4,2.0];
uniform = realSpaceInterlayerHopping(r, ...
    setInterlayerCompression(stack,0.05),1,2,1,2);
local = realSpaceInterlayerHopping(r,stack,1,2,1,2,[0.05,0.05,0.05]);
verifyEqual(testCase,local,uniform,'AbsTol',1e-15);
mixed = realSpaceInterlayerHopping(r,stack,1,2,1,2,[0,0.05,-0.1]);
verifyEqual(testCase,mixed(2),uniform(2),'AbsTol',1e-15);
verifyEqual(testCase,mixed(1), ...
    realSpaceInterlayerHopping(r(:,1),stack,1,2,1,2),'AbsTol',1e-15);
end

function h = legacyHopping(rVector,stack,j,k,alpha,beta)
% The constant-parameter implementation used before compression support.
rPhysical = sqrt(sum(rVector.^2,1));
r = rPhysical/stack.aG;
rHat = rVector./max(rPhysical,1e-14*stack.aG);
tauJHat = stack.tau(:,j)/norm(stack.tau(:,j));
tauKHat = stack.tau(:,k)/norm(stack.tau(:,k));
zJ = max(-1,min(1,sum(rHat.*tauJHat,1)));
zK = max(-1,min(1,sum(rHat.*tauKHat,1)));
T3J = 4*zJ.^3-3*zJ;
T3K = 4*zK.^3-3*zK;
T6J = 32*zJ.^6-48*zJ.^4+18*zJ.^2-1;
T6K = 32*zK.^6-48*zK.^4+18*zK.^2-1;
V0 = 0.3155*exp(-1.7543*r.^2).*cos(2.0010*r);
V3 = -0.0688*r.^2.*exp(-3.4692*(r-0.5212).^2);
V6 = -0.0083*exp(-2.8764*(r-1.5206).^2).*sin(1.5731*r);
h = V0+V3.*((-1)^(alpha+2)*T3J+(-1)^(beta+1)*T3K)+V6.*(T6J+T6K);
end
