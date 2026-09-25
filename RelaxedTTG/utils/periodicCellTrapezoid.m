function [b,weight] = periodicCellTrapezoid(A,N)
%PERIODICCELLTRAPEZOID Direct trapezoidal rule for a normalized cell average.
%
% avg_Gamma f(b) db = int_[0,1)^2 f(A*s) ds
%                     approximately (1/N^2) sum_mn f(A*[m;n]/N).
% This routine only returns quadrature nodes and weights; it uses no FFT.
assert(N >= 1 && N == fix(N),'N must be a positive integer.');

[s1,s2] = ndgrid((0:N-1)/N);
b = A*[s1(:),s2(:)].';
weight = ones(N^2,1)/N^2;
end
