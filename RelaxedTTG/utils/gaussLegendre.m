function [node,weight] = gaussLegendre(order)
%GAUSSLEGENDRE Nodes and weights on [-1,1] by Golub-Welsch.
assert(order >= 1 && order == fix(order), ...
    'order must be a positive integer.');

beta = 0.5./sqrt(1-(2*(1:order-1)).^-2);
T = diag(beta,1)+diag(beta,-1);
[V,D] = eig(T);
[node,index] = sort(diag(D));
weight = 2*V(1,index).^2;
node = node(:);
weight = weight(:);
end
