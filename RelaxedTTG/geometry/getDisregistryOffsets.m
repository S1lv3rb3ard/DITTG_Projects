function offsets = getDisregistryOffsets(stack,layer,tau)
%GETDISREGISTRYOFFSETS Configuration offsets D_j*tau of an orbital in layer j.
%
% offsets = getDisregistryOffsets(stack,layer,tau) returns a 2-by-3 array.
% Column t (t ~= layer) is (I - A_t*A_j^{-1})*tau and column t = layer is
% zero. An orbital at R+tau (R a lattice point of layer j) has the smooth
% disregistry configuration T_j*R + D_j*tau, so its relaxation displacement is
%
%   u_j( R + offsets(:,k), R + offsets(:,l) ),   [k,l] = setdiff(1:3,j,'stable').
%
% The raw configuration T_j*(R+tau) would shift the configuration by tau,
% i.e. by a sizeable fraction of a moire period in real space, and produces
% spurious bond strains of the order of the relaxation amplitude itself.
% See Zhu, Cazeaux, Luskin & Kaxiras, PRB 101, 224107 (2020), Eq. (5).
assert(ismember(layer,1:3),'layer must be 1, 2, or 3.');
assert(isequal(size(tau),[2,1]),'tau must be a 2-by-1 vector.');
offsets = zeros(2,3);
for t = setdiff(1:3,layer)
    offsets(:,t) = (eye(2)-stack.A{t}/stack.A{layer})*tau;
end
end
