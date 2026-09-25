function numDoF = getNumDoF(DoF)
%GETNUMDOF Number of reciprocal states retained in each represented layer.
layers = unique(DoF(:,9));
numDoF = zeros(numel(layers),1);
for k = 1:numel(layers)
    numDoF(k) = nnz(DoF(:,9) == layers(k));
end
end
