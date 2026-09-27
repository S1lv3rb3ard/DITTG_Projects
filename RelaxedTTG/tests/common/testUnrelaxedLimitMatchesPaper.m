function tests = testUnrelaxedLimitMatchesPaper
%TESTUNRELAXEDLIMITMATCHESPAPER Zero relaxation reproduces arXiv:2606.13434.
%
% With zero relaxation fields the cached relaxed Hamiltonian must equal the
% unrelaxed momentum-space Hamiltonian of Beard & Massatt (arXiv:2606.13434),
% assembled here directly from its formulas:
%
%   intralayer  (G,j alpha),(G,j beta):
%       sum_R exp(-i Q.(R + tau_alpha - tau_beta)) h_{alpha beta}(R + tau_alpha - tau_beta),
%       Q = q + G_k + G_l  (diagonal in G);
%   interlayer  (G',j alpha),(G'',k beta), nonzero only if G'_l = G''_l:
%       |Gamma_j Gamma_k|^(-1/2) exp(i G''_j.tau_{j alpha}) exp(-i G'_k.tau_{k beta})
%       * hhat_{j alpha,k beta}(q + G'_k + G'_l + G''_j),
%       hhat(xi) = int h(x) exp(-i xi.x) dx  (Fang-Kaxiras 2016 hopping).
%
% The continuous transform uses the same polar quadrature as the code, so the
% test checks the block placement, the Bloch phases and the hopping model.
tests = functiontests(localfunctions);
end

function testZeroRelaxationEqualsPaperHamiltonian(testCase)
a = 1.42*sqrt(3);
stack = getStack(a,[-1.4,0,2.8]);
DoF = getDoF(stack,0.35,4,'clean');
tA = [0,0.3302,0.23206,0.04969,-0.02499, ...
      0.00285,0.00204,-0.00014,-0.00029].';
tB = [-2.99251,-0.28983,0.02791,-0.00877,-0.01870, ...
      0.00621,-0.00256,-0.00018,-0.00033,-0.00264].';
shells = getInterpolatedIntralayerHoppingValues(stack,tA,tB,[]);
zeroField = @(xk,xl) zeros(size(xk),'like',xk);
fields = {zeroField,zeroField,zeroField};
q = stack.K(:,2)+[0.004;-0.003];

radialOrder = 40;
angularOrder = 64;
cutoff = 7*a;
options.verbose = false;
options.intralayer = struct('configurationGridSize',5,'useGPU',false, ...
    'verbose',false);
options.interlayer = struct( ...
    'realSpaceRadialOrder',radialOrder, ...
    'realSpaceAngularOrder',angularOrder, ...
    'realSpaceCutoff',cutoff, ...
    'configurationGridSize',5, ...
    'minimumConfigurationGridSize',3, ...
    'momentumInnerRadius',1e3, ...      % disable the smooth momentum cutoff
    'momentumOuterRadius',2e3, ...
    'useGPU',false, ...
    'verbose',false);
cache = prepareRelaxedHamiltonian(stack,DoF,shells,fields,q,options);
H = full(evaluateRelaxedHamiltonian(cache,q));

Href = paperHamiltonian(stack,DoF,shells,q,radialOrder,angularOrder,cutoff);
scale = max(abs(Href),[],'all');
verifyLessThan(testCase,max(abs(H-Href),[],'all')/scale,1e-10);
end

function H = paperHamiltonian(stack,DoF,shells,q,radialOrder,angularOrder,cutoff)
n = size(DoF,1);
H = complex(zeros(2*n));
layerOf = DoF(:,9);
G = @(row,layer) DoF(row,2*layer+(1:2)).';

% Intralayer blocks.
for i = 1:n
    j = layerOf(i);
    others = setdiff(1:3,j,'stable');
    Q = q+G(i,others(1))+G(i,others(2));
    shellMap = stack.A{j}/stack.A{2};
    RA = shellMap*shells.matA;           % lattice vectors R
    RB = shellMap*shells.matB;           % R + tau_B
    % Bonds r_row - r_col = R + tau_row - tau_col over the tabulated shells:
    % AA and BB: R; AB: -(R + tau_B) = R' + tau_A - tau_B; BA: R + tau_B.
    hAA = channelSum(Q,RA,shells.intraAA);
    hAB = channelSum(Q,-RB,shells.intraAB);
    hBA = channelSum(Q,RB,shells.intraAB);
    idx = 2*i-1:2*i;
    H(idx,idx) = [hAA,hAB;hBA,hAA];
end

% Interlayer blocks for the adjacent pairs (1,2) and (2,3).
[x,w] = polarQuadrature(radialOrder,angularOrder,cutoff);
normalization = @(j,k) sqrt(abs(det(stack.B{j}))*abs(det(stack.B{k})))/(2*pi)^2;
for pair = [1,2;2,3].'
    j = pair(1);
    k = pair(2);
    l = setdiff(1:3,[j,k]);
    rows = find(layerOf == j).';
    cols = find(layerOf == k).';
    for alpha = 1:2
        for beta = 1:2
            hx = realSpaceInterlayerHopping(x,stack,j,k,alpha,beta);
            tauJ = (alpha-1)*stack.tau(:,j);
            tauK = (beta-1)*stack.tau(:,k);
            for r = rows
                for c = cols
                    if norm(G(r,l)-G(c,l)) > 1e-9
                        continue
                    end
                    xi = q+G(r,k)+G(r,l)+G(c,j);
                    hhat = sum(w(:).'.*hx.*exp(-1i*(xi.'*x)));
                    value = normalization(j,k)*exp(1i*G(c,j).'*tauJ) ...
                        *exp(-1i*G(r,k).'*tauK)*hhat;
                    H(2*r-2+alpha,2*c-2+beta) = value;
                    H(2*c-2+beta,2*r-2+alpha) = conj(value);
                end
            end
        end
    end
end
end

function value = channelSum(Q,bonds,interpolant)
value = sum(exp(-1i*(Q.'*bonds)).*interpolant(vecnorm(bonds,2,1)));
end
