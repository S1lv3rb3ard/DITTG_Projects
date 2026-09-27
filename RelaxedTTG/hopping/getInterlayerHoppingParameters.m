function parameters = getInterlayerHoppingParameters(model,epsilon)
%GETINTERLAYERHOPPINGPARAMETERS Graphene interlayer p_z hopping parameters.
%
% parameters = getInterlayerHoppingParameters(model,epsilon) returns the ten
% parameters of the Fang-Kaxiras interlayer function
%
%   V0 = lambda0 exp(-xi0 r^2) cos(kappa0 r)
%   V3 = lambda3 r^2 exp(-xi3 (r-x3)^2)
%   V6 = lambda6 exp(-xi6 (r-x6)^2) sin(kappa6 r),   r = |r_parallel|/a,
%
% at the relative interlayer distance epsilon = d/d0 - 1 (negative under
% compression, d0 = 3.35 Angstrom). epsilon may be an array (for example
% one local distance per sample); every field then has the size of epsilon.
%
% Models
%   'carr2018'        (default) For -0.20 <= eps <= 0 the quadratic fits
%                     y(eps) = c0 + c1 eps + c2 eps^2 of Carr, Fang,
%                     Jarillo-Herrero & Kaxiras, PRB 98, 085144 (2018),
%                     Table I. For 0 < eps <= 0.15 (expansion, e.g. the AA
%                     regions of a corrugated moire) a C^1 continuation:
%                     the amplitudes lambda0, lambda3, lambda6 continue
%                     exponentially, lambda(eps) = lambda(0) exp(eps c1/c0),
%                     and the seven shape parameters linearly,
%                     y(eps) = c0 + c1 eps. The quadratic fits are not used for
%                     eps > 0 because they turn upward (lambda0 would have
%                     a minimum near eps = 0.12). The three amplitude decay
%                     rates at eps = 0 agree (-6.1, -5.9, -5.8), and match
%                     half the decay rate of the refined Kolmogorov-Crespi
%                     registry energy (-10.7), which scales as the square
%                     of the hopping. Outside [-0.20, 0.15] a warning is issued.
%   'fangKaxiras2016' the equilibrium parameters of Fang & Kaxiras, PRB 93,
%                     235153 (2016), used by RelaxedTTG before compression was
%                     added. Only epsilon = 0 is allowed.
%
% Carr et al. state eps = 1 - d/d0, but their pressure law
% P = A(exp(-B eps) - 1) and their fitted coefficients use eps < 0 for
% compression (eps = -0.1 gives 9.2 GPa and a larger lambda0); that is the
% convention used here.
if nargin < 1 || isempty(model)
    model = 'carr2018';
end
if nargin < 2 || isempty(epsilon)
    epsilon = 0;
end
assert(isnumeric(epsilon) && all(isfinite(epsilon(:))), ...
    'epsilon must be finite.');
names = {'lambda0','xi0','kappa0','lambda3','xi3','x3', ...
    'lambda6','xi6','x6','kappa6'};

switch lower(model)
    case 'carr2018'
        % Columns: c0, c1, c2 (eV for the lambda parameters).
        table = [ ...
             0.310, -1.882,  7.741;   % lambda0
             1.750, -1.618,  1.848;   % xi0
             1.990,  1.007,  2.427;   % kappa0
            -0.068,  0.399, -1.739;   % lambda3
             3.286, -0.914, 12.011;   % xi3
             0.500,  0.322,  0.908;   % x3
            -0.008,  0.046, -0.183;   % lambda6
             2.272, -0.721, -4.414;   % xi6
             1.217,  0.027, -0.658;   % x6
             1.562, -0.371, -0.134];  % kappa6
        amplitude = [true,false,false,true,false,false,true,false,false,false];
        if any(epsilon(:) < -0.20-1e-12) || any(epsilon(:) > 0.15+1e-12)
            warning('getInterlayerHoppingParameters:extrapolation', ...
                ['epsilon in [%.4g, %.4g] extends beyond [-0.20, 0.15]: ' ...
                 'the fit (eps <= 0) and its continuation (eps > 0) ' ...
                 'are extrapolated.'],min(epsilon(:)),max(epsilon(:)));
        end
        compressed = epsilon <= 0;
        values = cell(1,10);
        for index = 1:10
            c = table(index,:);
            value = c(1)+c(2)*epsilon+c(3)*epsilon.^2;
            if amplitude(index)
                continued = c(1)*exp(epsilon*c(2)/c(1));
            else
                continued = c(1)+c(2)*epsilon;
            end
            value(~compressed) = continued(~compressed);
            values{index} = value;
        end
    case {'fangkaxiras2016','fang-kaxiras2016'}
        assert(all(abs(epsilon(:)) < 1e-14), ...
            ['The Fang-Kaxiras (2016) parameters have no compression ' ...
             'dependence. Use the ''carr2018'' model for epsilon ~= 0.']);
        constants = [0.3155,1.7543,2.0010,-0.0688,3.4692,0.5212, ...
            -0.0083,2.8764,1.5206,1.5731];
        values = arrayfun(@(v) v+zeros(size(epsilon)),constants, ...
            'UniformOutput',false);
    otherwise
        error('Unknown interlayer hopping model "%s".',model);
end

parameters = cell2struct(values(:),names(:),1);
parameters.model = lower(model);
parameters.epsilon = epsilon;
end
