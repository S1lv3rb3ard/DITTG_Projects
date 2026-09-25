function basis = getMoireBasis(stack,DoF,W)
%GETMOIREBASIS Build layer-resolved projectors for total-DoS quadrature.
invB = cellfun(@inv,stack.B,'UniformOutput',false);
integerComponents = round(blkdiag(invB{:})*DoF(:,3:8).');
physicalG = cat(2,stack.B{:})*integerComponents;

basis = cell(1,3);
for layer = 1:3
    if layer == 1
        C = stack.B{1}-stack.B{2};
    elseif layer == 2
        basis{2} = basis{1};
        continue
    else
        C = stack.B{2}-stack.B{3};
    end
    limitN = floor(W*vecnorm(inv(C),2,2));
    [n1,n2] = ndgrid(-limitN(1):limitN(1),-limitN(2):limitN(2));
    validG = C*[n1(:),n2(:)].';
    validG = validG(:,vecnorm(validG,2,1) < W);

    % Avoid a Statistics Toolbox dependency on pdist2.
    distanceSquared = reshape(sum(validG.^2,1),[],1) ...
        +sum(physicalG.^2,1)-2*validG.'*physicalG;
    [~,column] = find(abs(distanceSquared) < 1e-28);
    orbital = reshape([2*(column-1)+1;2*(column-1)+2],[],1);
    basis{layer} = sparse(orbital,1:numel(orbital),1, ...
        2*size(DoF,1),numel(orbital));
end
end
