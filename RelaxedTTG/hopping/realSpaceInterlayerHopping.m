function h = realSpaceInterlayerHopping(rVector,stack,j,k,alpha,beta,localEpsilon)
%REALSPACEINTERLAYERHOPPING Appendix-A graphene hopping h_{k beta}^{j alpha}.
% rVector is a 2-by-N array of physical relative positions. The ten radial
% and angular parameters are those of the interface (j,k) at its uniform
% compression, stack.interlayer (see setInterlayerCompression); stacks
% without that field use the 'carr2018' parameters at zero compression.
%
% Optional localEpsilon (scalar or 1-by-N) gives the relative interlayer
% distance d/d0 - 1 of each sample, for example from out-of-plane
% corrugation; it replaces the uniform compression of the interface and
% uses the interface's hopping model.
rPhysical = sqrt(sum(rVector.^2,1));
r = rPhysical/stack.aG;
safeNorm = max(rPhysical,1e-14*stack.aG);
rHat = rVector./safeNorm;

tauJHat = stack.tau(:,j)/norm(stack.tau(:,j));
tauKHat = stack.tau(:,k)/norm(stack.tau(:,k));
tauJHat = toDevice(tauJHat,isa(rVector,'gpuArray'));
tauKHat = toDevice(tauKHat,isa(rVector,'gpuArray'));

zJ = max(-1,min(1,sum(rHat.*tauJHat,1)));
zK = max(-1,min(1,sum(rHat.*tauKHat,1)));
T3J = 4*zJ.^3-3*zJ;
T3K = 4*zK.^3-3*zK;
T6J = 32*zJ.^6-48*zJ.^4+18*zJ.^2-1;
T6K = 32*zK.^6-48*zK.^4+18*zK.^2-1;

if nargin < 7 || isempty(localEpsilon)
    p = interlayerParametersForPair(stack,j,k);
else
    assert(isscalar(localEpsilon) || isequal(size(localEpsilon),size(r)), ...
        'localEpsilon must be a scalar or match the 1-by-N samples.');
    p = getInterlayerHoppingParameters(interlayerModel(stack),localEpsilon);
end
% Elementwise throughout: the parameters are scalars or 1-by-N arrays.
V0 = p.lambda0.*exp(-p.xi0.*r.^2).*cos(p.kappa0.*r);
V3 = p.lambda3.*r.^2.*exp(-p.xi3.*(r-p.x3).^2);
V6 = p.lambda6.*exp(-p.xi6.*(r-p.x6).^2).*sin(p.kappa6.*r);

theta3 = (-1)^(alpha+2)*T3J+(-1)^(beta+1)*T3K;
theta6 = T6J+T6K;
h = V0+V3.*theta3+V6.*theta6;
end

function p = interlayerParametersForPair(stack,j,k)
pair = sort([j,k]);
if isequal(pair,[1,2])
    index = 1;
elseif isequal(pair,[2,3])
    index = 2;
else
    error('realSpaceInterlayerHopping:unsupportedPair', ...
        ['Only the adjacent interfaces (1,2) and (2,3) are ' ...
        'parameterized; the pair (%d,%d) is not.'],j,k);
end
if isfield(stack,'interlayer') && isfield(stack.interlayer,'parameters')
    p = stack.interlayer.parameters{index};
else
    p = getInterlayerHoppingParameters('carr2018',0);
end
end

function model = interlayerModel(stack)
if isfield(stack,'interlayer') && isfield(stack.interlayer,'model')
    model = stack.interlayer.model;
else
    model = 'carr2018';
end
end
