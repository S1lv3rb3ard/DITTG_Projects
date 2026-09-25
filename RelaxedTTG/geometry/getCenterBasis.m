function X = getCenterBasis(DoF,x)
%GETCENTERBASIS Select the (G=0,j,alpha) momentum-space orbitals.
%
% The six center orbitals are ordered by the DoF ordering and then by the
% two sublattices.  Selecting x=1:6 produces the projector used in the
% averaged momentum LDoS of equations (3.20) and (4.51) of the paper.
centerState = find(all(DoF(:,3:8) == 0,2)).';
centerOrbital = reshape([2*centerState-1;2*centerState],[],1);
assert(numel(centerOrbital) == 6, ...
    'Expected exactly two center orbitals in each of three layers.');
if nargin < 2 || isempty(x)
    x = 1:6;
end
assert(all(x >= 1 & x <= 6 & x == fix(x)), ...
    'x must contain indices from 1 through 6.');
X = sparse(centerOrbital(x),1:numel(x),1, ...
    2*size(DoF,1),numel(x));
end
