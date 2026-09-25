function T = Chebyshev(P,scaledE,mode)
%CHEBYSHEV Evaluate T_0 through T_P at the scaled energies.
assert(P >= 1 && P == fix(P),'P must be a positive integer.');
assert(all(abs(scaledE) < 1,'all'), ...
    'Every scaled energy must lie strictly inside (-1,1).');
T = zeros(P+1,numel(scaledE),'like',scaledE);
if strcmpi(mode,'gpu')
    T = gpuArray(T);
    scaledE = gpuArray(scaledE);
end
T(1,:) = 1;
T(2,:) = scaledE;
twiceEnergy = 2.*scaledE;
for order = 3:P+1
    T(order,:) = twiceEnergy.*T(order-1,:)-T(order-2,:);
end
end
