% Dumps the zero-relaxation RelaxedTTG Hamiltonian (paper model, arXiv:2606.13434)
% at a few momenta, for comparison with an independent Python port.
root = fileparts(mfilename('fullpath'));
addpath(genpath(root));
a = 1.42*sqrt(3);
stack = getStack(a,[-1.5,0,2.0]);
W = 0.35; L = 10;
DoF = getDoF(stack,W,L,'clean');
tA = [0,0.3302,0.23206,0.04969,-0.02499,0.00285,0.00204,-0.00014,-0.00029].';
tB = [-2.99251,-0.28983,0.02791,-0.00877,-0.01870,0.00621,-0.00256,-0.00018,-0.00033,-0.00264].';
shells = getInterpolatedIntralayerHoppingValues(stack,tA,tB,[]);
zeroField = @(xk,xl) zeros(size(xk),'like',xk);
fields = {zeroField,zeroField,zeroField};
qs = [stack.K(:,1)+[0.03;-0.02], stack.K(:,2)+[-0.02;0.01], stack.K(:,3)+[0.015;0.03]];
opts.verbose = false;
opts.intralayer = struct('configurationGridSize',5,'useGPU',false,'verbose',false);
opts.interlayer = struct('realSpaceRadialOrder',80,'realSpaceAngularOrder',128, ...
    'realSpaceCutoff',8*a,'configurationGridSize',5,'minimumConfigurationGridSize',3, ...
    'momentumInnerRadius',1e3,'momentumOuterRadius',2e3,'useGPU',false,'verbose',false);
optsDefault = opts;
optsDefault.interlayer = rmfield(optsDefault.interlayer,{'momentumInnerRadius','momentumOuterRadius'});
H = cell(1,3); Hdefault = cell(1,3);
for i = 1:3
    cache = prepareRelaxedHamiltonian(stack,DoF,shells,fields,qs(:,i),opts);
    H{i} = evaluateRelaxedHamiltonian(cache,qs(:,i));
    cache = prepareRelaxedHamiltonian(stack,DoF,shells,fields,qs(:,i),optsDefault);
    Hdefault{i} = evaluateRelaxedHamiltonian(cache,qs(:,i));
end
H1 = H{1}; H2 = H{2}; H3 = H{3}; D1 = Hdefault{1}; D2 = Hdefault{2}; D3 = Hdefault{3};
save(fullfile(root,'validate_tb_port.mat'),'DoF','qs','H1','H2','H3','D1','D2','D3','-v7');
fprintf('validate_tb_port: saved %d x %d Hamiltonians\n',size(H1,1),size(H1,2));
