function h = realSpaceInterlayerHopping(rVector,stack,j,k,alpha,beta)
%REALSPACEINTERLAYERHOPPING Appendix-A graphene hopping h_{k beta}^{j alpha}.
% rVector is a 2-by-N array of physical relative positions.
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

V0 = 0.3155*exp(-1.7543*r.^2).*cos(2.0010*r);
V3 = -0.0688*r.^2.*exp(-3.4692*(r-0.5212).^2);
V6 = -0.0083*exp(-2.8764*(r-1.5206).^2).*sin(1.5731*r);

theta3 = (-1)^(alpha+2)*T3J+(-1)^(beta+1)*T3K;
theta6 = T6J+T6K;
h = V0+V3.*theta3+V6.*theta6;
end
