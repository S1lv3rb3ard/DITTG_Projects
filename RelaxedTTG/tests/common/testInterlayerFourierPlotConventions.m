function tests = testInterlayerFourierPlotConventions
%TESTINTERLAYERFOURIERPLOTCONVENTIONS Verify slice/mode and normalization.
tests = functiontests(localfunctions);
end

function spectatorIndependentLimit(testCase)
stack = getStack(1.42*sqrt(3),[-0.3,0,0.3]);
zeroFields = cell(1,3);
for layer = 1:3
    zeroFields{layer} = @(xk,xl) zeros(size(xk),'like',xk); %#ok<INUSD>
end
Q = [0,0.5,-0.4;0,0.2,0.7];
baseOptions = struct( ...
    'realSpaceRadialOrder',12, ...
    'realSpaceAngularOrder',18, ...
    'realSpaceCutoff',8, ...
    'configurationGridSize',5, ...
    'minimumConfigurationGridSize',5, ...
    'batchSize',3, ...
    'useGPU',false, ...
    'spectatorFraction',[0.23;0.41], ...
    'includeBlochNormalization',false);

sliceOptions = baseOptions;
sliceOptions.transformType = 'spectatorSlice';
sliceValue = computeRelaxedInterlayerFourierHopping( ...
    stack,zeroFields,[1,2],1,1,Q,[0;0],sliceOptions);

modeOptions = baseOptions;
modeOptions.transformType = 'spectatorMode';
modeZero = computeRelaxedInterlayerFourierHopping( ...
    stack,zeroFields,[1,2],1,1,Q,[0;0],modeOptions);
modeNonzero = computeRelaxedInterlayerFourierHopping( ...
    stack,zeroFields,[1,2],1,1,Q,[1;0],modeOptions);

verifyLessThan(testCase,max(abs(sliceValue-modeZero)),1e-12);
verifyLessThan(testCase,max(abs(modeNonzero)),1e-12);

normalizedOptions = sliceOptions;
normalizedOptions.includeBlochNormalization = true;
normalized = computeRelaxedInterlayerFourierHopping( ...
    stack,zeroFields,[1,2],1,1,Q,[0;0],normalizedOptions);
factor = sqrt(abs(det(stack.B{1}))*abs(det(stack.B{2})));
verifyLessThan(testCase,max(abs(normalized-factor*sliceValue)),1e-12);
end
