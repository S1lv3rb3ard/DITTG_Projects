function [relaxationFields,info] = loadJuliaTrilayerRelaxation( ...
    minimizerFile,dataFile,stack,options)
%LOADJULIATRILAYERRELAXATION Relaxation fields from the Julia minimizer.
%
% [relaxationFields,info] = loadJuliaTrilayerRelaxation(minimizerFile, ...
%     dataFile,stack,options)
%
% Reads the raw binary files written by example2.jl,
%   triG_data_<t12>_<t23>_<N>.jld       N, theta1, theta3, E, P, K, G
%   triG_minimizer_<t12>_<t23>_<N>.jld  u, a 2-by-N-by-N-by-N-by-N-by-3 array
% and returns relaxationFields{j}(xk,xl) in the convention used throughout
% RelaxedTTG: [k,l] = setdiff(1:3,j,'stable'), xk and xl are 2-by-M physical
% coordinates of a point of layer j modulo the lattices of layers k and l,
% and the output is the 2-by-M in-plane displacement in Angstrom.
%
% Conventions reconciled here (see README, "Relaxation fields"):
%  * Julia stores u on the hull in fractional coordinates of the shift
%    b = R^(t) - x (Zhu et al., PRB 101, 224107), i.e. MINUS the position
%    modulo the lattice. RelaxedTTG uses +x modulo the lattice.
%  * Julia layer j stores the configuration pairs relative to layers
%    (2,3), (1,3), (2,1) for j = 1,2,3. Layer 3 is therefore reordered.
%  * Julia rotates layers by [cos t, sin t; -sin t, cos t] and uses a
%    lattice basis at +-30 degrees; getStack uses counterclockwise rotations
%    and a basis at 0/60 degrees. A frame rotation Phi (30 degrees plus the
%    stack angle of layer 2) maps Julia coordinates and vectors to RelaxedTTG.
%    It is found automatically and checked against every layer and orbital.
%
% The hull data are converted to a trigonometric (spectral) interpolant,
%   u_j(xk,xl) = Re sum_m c_m exp(i(Gk_m.xk + Gl_m.xl)),
% with Gk_m, Gl_m in the reciprocal lattices of layers k and l, so the
% fields are exactly periodic and reproduce the Julia grid values at the
% grid configurations. Modes with amplitude below
% options.relativeTolerance*max are discarded.
%
% Options: relativeTolerance (1e-6), maximumModes (inf), chunkElements
% (2e7), angleTolerance (1e-8 in fractional lattice coordinates), verbose.
if nargin < 4
    options = struct;
end
defaults = struct( ...
    'relativeTolerance',1e-6, ...
    'maximumModes',inf, ...
    'chunkElements',2e7, ...
    'latticeTolerance',1e-7, ...
    'verbose',true);
options = mergeOptions(options,defaults);

julia = readJuliaRelaxationFiles(minimizerFile,dataFile);
N = julia.N;
jrot = @(t) [cos(t),sin(t);-sin(t),cos(t)];
thetaJulia = [julia.theta1,0,julia.theta3];
tE = cell(1,3);
tP = cell(1,3);
for layer = 1:3
    tE{layer} = jrot(thetaJulia(layer))*julia.E;
    tP{layer} = jrot(thetaJulia(layer))*julia.P;
end

Phi = findFrameRotation(stack,tE,tP,options.latticeTolerance);
if isempty(Phi)
    expected = stack.theta(2)-rad2deg([julia.theta1,0,julia.theta3]);
    error('loadJuliaTrilayerRelaxation:geometryMismatch', ...
        ['The stack does not match the Julia geometry (theta12 = %.6g, ' ...
         'theta23 = %.6g degrees, |a| = %.6g). Use getStack(%.6g,[%s]) ' ...
         'or the same angles shifted by a common rotation.'], ...
        rad2deg(julia.theta1),rad2deg(julia.theta3), ...
        norm(julia.E(:,1)),norm(julia.E(:,1)),num2str(expected,'%.6g '));
end

% Julia configuration pairs for each layer: dims (1,2) and (3,4) of the hull.
juliaPairs = {[2,3],[1,3],[2,1]};
frequency = [0:ceil(N/2)-1,-floor(N/2):-1];
[n1,n2,n3,n4] = ndgrid(frequency);
allModes = [n1(:),n2(:),n3(:),n4(:)];
clear n1 n2 n3 n4

relaxationFields = cell(1,3);
info.N = N;
info.thetaJuliaDegrees = rad2deg([julia.theta1,julia.theta3]);
info.Phi = Phi;
info.frameAngleDegrees = atan2d(Phi(2,1),Phi(1,1));
info.bulkModulus = julia.K;
info.shearModulus = julia.G;
info.numModes = zeros(1,3);
info.maxDisplacement = zeros(1,3);
info.discardedAmplitude = zeros(1,3);
for layer = 1:3
    otherLayers = setdiff(1:3,layer,'stable');
    Ux = fftn(reshape(julia.u(1,:,:,:,:,layer),[N,N,N,N]))/N^4;
    Uy = fftn(reshape(julia.u(2,:,:,:,:,layer),[N,N,N,N]))/N^4;
    coefficients = [Ux(:),Uy(:)];
    modes = allModes;
    magnitude = max(abs(coefficients),[],2);
    keep = magnitude >= options.relativeTolerance*max(magnitude);
    if isfinite(options.maximumModes) && nnz(keep) > options.maximumModes
        [~,order] = sort(magnitude,'descend');
        keep = false(size(keep));
        keep(order(1:options.maximumModes)) = true;
    end
    info.discardedAmplitude(layer) = sum(sum(abs(coefficients(~keep,:))));
    coefficients = coefficients(keep,:);
    modes = modes(keep,:);
    [modes,coefficients] = splitNyquistModes(modes,coefficients,N);

    % Physical wave vectors in RelaxedTTG coordinates. The Julia phase is
    % 2*pi*n.gamma with gamma = tE_t^{-1}(-x_J) and x_J = Phi^{-1} x.
    wave = cell(1,2);
    for slot = 1:2
        relative = otherLayers(slot);
        pairIndex = find(juliaPairs{layer} == relative);
        assert(numel(pairIndex) == 1,'Internal pair ordering error.');
        n = modes(:,2*pairIndex-1:2*pairIndex).';
        wave{slot} = -Phi*(2*pi*inv(tE{relative}).')*n;
    end
    field.layer = layer;
    field.Gk = wave{1};
    field.Gl = wave{2};
    field.coefficients = (Phi*coefficients.').';
    field.chunkElements = options.chunkElements;
    relaxationFields{layer} = @(xk,xl) evaluateFourierRelaxationField( ...
        field,xk,xl);
    info.numModes(layer) = size(modes,1);
    info.maxDisplacement(layer) = max(abs(julia.u(:,:,:,:,:,layer)),[],'all');
    info.field{layer} = field;
end

if options.verbose
    fprintf(['Loaded Julia relaxation: theta12 = %.4g, theta23 = %.4g deg, ' ...
        'N = %d, frame rotation %.4g deg.\n'], ...
        info.thetaJuliaDegrees(1),info.thetaJuliaDegrees(2),N, ...
        info.frameAngleDegrees);
    for layer = 1:3
        fprintf(['  layer %d: %d Fourier modes, max |u| = %.4g A, ' ...
            'discarded amplitude <= %.2g A\n'],layer, ...
            info.numModes(layer),info.maxDisplacement(layer), ...
            info.discardedAmplitude(layer));
    end
end
end

function Phi = findFrameRotation(stack,tE,tP,tolerance)
Phi = [];
for angle = stack.theta(2)+(0:30:330)
    candidate = rotationMatrixDegrees(angle);
    ok = true;
    for layer = 1:3
        M = stack.A{layer}\(candidate*tE{layer});
        shift = stack.A{layer}\(candidate*tP{layer}-stack.tau(:,layer));
        ok = ok && all(abs(M-round(M)) < tolerance,'all') && ...
            abs(abs(det(round(M)))-1) < 0.5 && ...
            all(abs(shift-round(shift)) < tolerance);
    end
    if ok
        Phi = candidate;
        return
    end
end
end

function [modes,coefficients] = splitNyquistModes(modes,coefficients,N)
% For even N the Nyquist frequency -N/2 is shared equally with +N/2 so the
% interpolant is real and reproduces the grid values.
if mod(N,2) ~= 0
    return
end
for dim = 1:4
    nyquist = modes(:,dim) == -N/2;
    if any(nyquist)
        extraModes = modes(nyquist,:);
        extraModes(:,dim) = N/2;
        coefficients(nyquist,:) = 0.5*coefficients(nyquist,:);
        modes = [modes;extraModes]; %#ok<AGROW>
        coefficients = [coefficients;coefficients(nyquist,:)]; %#ok<AGROW>
    end
end
end
