function stack = setInterlayerCompression(stack,epsilon,model)
%SETINTERLAYERCOMPRESSION Uniform vertical compression of each interface.
%
% stack = setInterlayerCompression(stack,epsilon,model)
%
% epsilon = [eps12,eps23] (or a scalar for both interfaces) is the relative
% interlayer distance eps = d/d0 - 1, negative under compression, with
% d0 = 3.35 Angstrom. model selects the interlayer hopping parameterization
% (see getInterlayerHoppingParameters; default 'carr2018'). The interlayer
% hopping of the adjacent pairs (1,2) and (2,3) is evaluated with these
% parameters wherever realSpaceInterlayerHopping is called.
%
% To prescribe a uniaxial pressure P (GPa) instead, use
%   stack = setInterlayerCompression(stack,compressionFromPressure(P));
if nargin < 2 || isempty(epsilon)
    epsilon = 0;
end
if nargin < 3 || isempty(model)
    model = 'carr2018';
end
if isscalar(epsilon)
    epsilon = [epsilon,epsilon];
end
assert(numel(epsilon) == 2 && all(isfinite(epsilon)), ...
    'epsilon must be a scalar or a two-element vector [eps12,eps23].');
assert(all(epsilon > -1),'epsilon must exceed -1.');

d0 = 3.35;
stack.interlayer.model = lower(model);
stack.interlayer.d0 = d0;
stack.interlayer.compression = epsilon(:).';
stack.interlayer.spacing = d0*(1+epsilon(:).');
stack.interlayer.parameters = { ...
    getInterlayerHoppingParameters(model,epsilon(1)), ...
    getInterlayerHoppingParameters(model,epsilon(2))};
end
