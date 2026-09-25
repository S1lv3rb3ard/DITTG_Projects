function shells = getInterpolatedIntralayerHoppingValues(stack,tA,tB,rCut)
%GETINTERPOLATEDINTRALAYERHOPPINGVALUES Build shell hopping interpolants.
a = stack.aG;
mMax = max(20,2*max(numel(tA),numel(tB)));
[m1,m2] = ndgrid(0:mMax);
vA = unique(m1(:).^2+m1(:).*m2(:)+m2(:).^2,'sorted');
vB = unique(3*(m1(:).^2+m1(:).*m2(:)+m2(:).^2) ...
    +3*(m1(:)+m2(:))+1,'sorted');

rA = a*sqrt(vA(1:numel(tA)));
rB = (a/sqrt(3))*sqrt(vB(1:numel(tB)));
shells.rA = rA;
shells.rB = rB;

if nargin < 4 || isempty(rCut)
    rCutA = rA(end)+(rA(end)-rA(end-1));
    rCutB = rB(end)+(rB(end)-rB(end-1));
elseif isscalar(rCut)
    rCutA = rCut;
    rCutB = rCut;
else
    rCutA = rCut(1);
    rCutB = rCut(2);
end
assert(rCutA > rA(end),'The A cutoff must exceed the final A-shell radius.');
assert(rCutB > rB(end),'The B cutoff must exceed the final B-shell radius.');

shells.intraAA = griddedInterpolant( ...
    [rA;rCutA],[tA(:);0],'pchip','linear');
shells.intraAB = griddedInterpolant( ...
    [rB;rCutB],[tB(:);0],'pchip','linear');

[x1,x2] = ndgrid(-10:10);
integerSites = [x1(:),x2(:)].';
sitesA = stack.A{2}*integerSites;
sitesB = sitesA+stack.tau(:,2);
distanceA = vecnorm(sitesA,2,1);
distanceB = vecnorm(sitesB,2,1);
tolerance = 1e-10*a;

shells.A = arrayfun(@(r) sitesA(:,abs(distanceA-r)<tolerance), ...
    rA,'UniformOutput',false);
shells.B = arrayfun(@(r) sitesB(:,abs(distanceB-r)<tolerance), ...
    rB,'UniformOutput',false);
shells.matA = cat(2,shells.A{:});
shells.matB = cat(2,shells.B{:});
end
