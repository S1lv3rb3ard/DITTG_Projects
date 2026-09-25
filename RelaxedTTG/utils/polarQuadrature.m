function [x,weight] = polarQuadrature(radialOrder,angularOrder,cutoff)
%POLARQUADRATURE Direct quadrature for a disk-truncated R^2 integral.
% The returned weights include the polar Jacobian r.
[radialNode,radialWeight] = gaussLegendre(radialOrder);
radius = 0.5*cutoff*(radialNode+1);
radialWeight = 0.5*cutoff*radialWeight.*radius;

[angularNode,angularWeight] = gaussLegendre(angularOrder);
angle = pi*(angularNode+1);
angularWeight = pi*angularWeight;

[radiusGrid,angleGrid] = ndgrid(radius,angle);
[radialWeightGrid,angularWeightGrid] = ndgrid( ...
    radialWeight,angularWeight);
x = [radiusGrid(:).*cos(angleGrid(:)), ...
     radiusGrid(:).*sin(angleGrid(:))].';
weight = radialWeightGrid(:).*angularWeightGrid(:);
end
