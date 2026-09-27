function channels = sampleRelaxedIntralayerChannels( ...
    stack,shellsReference,relaxation,layer,bk,bl,options)
%SAMPLERELAXEDINTRALAYERCHANNELS Relaxed intralayer hoppings on hull samples.
%
% channels = sampleRelaxedIntralayerChannels(stack,shellsReference, ...
%     relaxation,layer,bk,bl,options)
%
% bk and bl are 2-by-Nc configuration samples of the COLUMN lattice point,
% b = T_j R' modulo the lattices of the complementary layers
% [k,l] = setdiff(1:3,layer,'stable'). The four channels are returned in the
% order AA, AB, BA, BB, where the first letter is the row orbital and the
% second the column orbital. For each channel:
%
%   bond0        2-by-S unrelaxed bonds  r_row - r_col = d + tau_row - tau_col
%   latticeShift 2-by-S lattice differences d = R - R'
%   tauRow/tauCol orbital offsets
%   relaxedLength Nc-by-S |bond0 + u_row(T_j d + b) - u_col(b)|
%   hopping      Nc-by-S hopping values at the relaxed lengths
%   relaxedBond  2-by-S relaxed bond vectors (only when Nc = 1)
%
% Orbital fields are evaluated at the disregistry-consistent configuration
% u_{j,alpha}(b) = u_j(b + D_j tau_alpha); see getDisregistryOffsets.
if nargin < 7
    options = struct;
end
defaults = struct('shellReferenceLayer',2,'useGPU',false);
options = mergeOptions(options,defaults);
assert(ismember(layer,1:3),'layer must be 1, 2, or 3.');
assert(size(bk,1) == 2 && isequal(size(bk),size(bl)), ...
    'bk and bl must be matching 2-by-N arrays.');

otherLayers = setdiff(1:3,layer,'stable');
kLayer = otherLayers(1);
lLayer = otherLayers(2);
shellMap = stack.A{layer}/stack.A{options.shellReferenceLayer};
RA = shellMap*shellsReference.matA;     % lattice vectors d
RB = shellMap*shellsReference.matB;     % vectors d + tau
t = stack.tau(:,layer);
tauOrbital = [zeros(2,1),t];
offsetA = getDisregistryOffsets(stack,layer,tauOrbital(:,1));
offsetB = getDisregistryOffsets(stack,layer,tauOrbital(:,2));
offsets = {offsetA,offsetB};

names = {'AA','AB','BA','BB'};
rowOrbital = [1,1,2,2];
colOrbital = [1,2,1,2];
bonds = {RA,-RB,RB,RA};
interpolants = {shellsReference.intraAA,shellsReference.intraAB, ...
    shellsReference.intraAB,shellsReference.intraAA};

useGPU = options.useGPU;
numberConfig = size(bk,2);
channels = struct('name',names,'rowOrbital',num2cell(rowOrbital), ...
    'colOrbital',num2cell(colOrbital));
for c = 1:4
    a = rowOrbital(c);
    b = colOrbital(c);
    bond0 = bonds{c};
    latticeShift = bond0-tauOrbital(:,a)+tauOrbital(:,b);
    offRow = toDevice(offsets{a},useGPU);
    offCol = toDevice(offsets{b},useGPU);
    uCol = relaxation(bk+offCol(:,kLayer),bl+offCol(:,lLayer));
    numberShells = size(bond0,2);
    lengths = zeros(numberConfig,numberShells,'like',bk);
    relaxedBonds = zeros(2,numberShells*(numberConfig == 1),'like',bk);
    for s = 1:numberShells
        d = toDevice(latticeShift(:,s),useGPU);
        uRow = relaxation(bk+d+offRow(:,kLayer),bl+d+offRow(:,lLayer));
        relaxedBond = toDevice(bond0(:,s),useGPU)+uRow-uCol;
        lengths(:,s) = vecnorm(relaxedBond,2,1).';
        if numberConfig == 1
            relaxedBonds(:,s) = relaxedBond;
        end
    end
    channels(c).bond0 = bond0;
    channels(c).latticeShift = latticeShift;
    channels(c).tauRow = tauOrbital(:,a);
    channels(c).tauCol = tauOrbital(:,b);
    channels(c).relaxedLength = lengths;
    channels(c).relaxedBond = toHost(relaxedBonds,useGPU);
    channels(c).hopping = toDevice( ...
        interpolants{c}(toHost(lengths,useGPU)),useGPU);
end
end
