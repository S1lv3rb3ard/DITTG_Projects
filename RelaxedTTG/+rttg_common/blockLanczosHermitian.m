function lanczos = blockLanczosHermitian(H,V,options)
%BLOCKLANCZOSHERMITIAN Fully reorthogonalized Hermitian block Lanczos.
%
% lanczos = blockLanczosHermitian(H,V,options) constructs an orthonormal
% block-Krylov basis for H starting from the columns of V.  Every new block
% is orthogonalized against the complete accumulated basis.  Two passes are
% used by default; this is the orthogonal correction that suppresses ghost
% Ritz values in deep recursions.
%
% The returned projectedMatrix is explicitly Hermitian.  Consequently its
% compressed resolvent is a matrix-valued Herglotz function even when the
% recursion is stopped at finite depth.
if nargin < 3
    options = struct;
end
defaults = struct( ...
    'maximumBlockSteps',100, ...
    'reorthogonalizationPasses',2, ...
    'breakdownTolerance',1e-12, ...
    'hermiticityTolerance',1e-12, ...
    'verbose',false);
options = mergeOptions(options,defaults);

assert(ismatrix(H) && size(H,1) == size(H,2), ...
    'H must be a square matrix.');
assert(size(V,1) == size(H,1) && ~isempty(V), ...
    'V must be nonempty and have the same number of rows as H.');
assert(options.maximumBlockSteps >= 1 && ...
    options.maximumBlockSteps == fix(options.maximumBlockSteps), ...
    'maximumBlockSteps must be a positive integer.');
assert(options.reorthogonalizationPasses >= 1 && ...
    options.reorthogonalizationPasses == ...
    fix(options.reorthogonalizationPasses), ...
    'reorthogonalizationPasses must be a positive integer.');

matrixNorm = norm(H,'fro');
hermiticityError = norm(H-H','fro')/max(matrixNorm,eps);
assert(hermiticityError <= options.hermiticityTolerance, ...
    'H must be Hermitian to the requested tolerance.');

% Rank-revealing factorization V=Q0*R0 also supports non-orthonormal V.
[leftVector,singularValue,rightVector] = svd(full(V),'econ');
singularValue = diag(singularValue);
assert(~isempty(singularValue) && singularValue(1) > 0, ...
    'V must have at least one nonzero column.');
retained = singularValue > ...
    options.breakdownTolerance*singularValue(1);
Q0 = leftVector(:,retained);
R0 = diag(singularValue(retained))*rightVector(:,retained)';

blocks = {Q0};
basis = Q0;
blockRanges = {1:size(Q0,2)};
breakdown = false;

for step = 1:options.maximumBlockSteps-1
    residual = H*blocks{end};
    for pass = 1:options.reorthogonalizationPasses
        residual = residual-basis*(basis'*residual);
    end

    [nextVector,nextValue,~] = svd(full(residual),'econ');
    nextValue = diag(nextValue);
    if isempty(nextValue) || nextValue(1) == 0
        breakdown = true;
        break
    end
    retained = nextValue > ...
        options.breakdownTolerance*max(matrixNorm,nextValue(1));
    remainingDimension = size(H,1)-size(basis,2);
    retainedIndex = find(retained,remainingDimension,'first');
    if isempty(retainedIndex)
        breakdown = true;
        break
    end

    nextBlock = nextVector(:,retainedIndex);
    % Correct once more after rank truncation.
    for pass = 1:options.reorthogonalizationPasses
        nextBlock = nextBlock-basis*(basis'*nextBlock);
    end
    [nextBlock,~] = qr(nextBlock,0);

    first = size(basis,2)+1;
    basis = [basis,nextBlock]; %#ok<AGROW>
    blocks{end+1} = nextBlock; %#ok<AGROW>
    blockRanges{end+1} = first:size(basis,2); %#ok<AGROW>
    if size(basis,2) == size(H,1)
        breakdown = true;
        break
    end
end

rawProjectedMatrix = full(basis'*(H*basis));
rawProjectedMatrix = 0.5*(rawProjectedMatrix+rawProjectedMatrix');
blockJacobiMatrix = zeros(size(rawProjectedMatrix));
for rowBlock = 1:numel(blockRanges)
    for columnBlock = max(1,rowBlock-1):min(numel(blockRanges),rowBlock+1)
        rows = blockRanges{rowBlock};
        columns = blockRanges{columnBlock};
        blockJacobiMatrix(rows,columns) = ...
            rawProjectedMatrix(rows,columns);
    end
end
blockJacobiMatrix = 0.5*(blockJacobiMatrix+blockJacobiMatrix');
[ritzVector,ritzValue] = eig(blockJacobiMatrix,'vector');
[ritzValue,order] = sort(real(ritzValue),'ascend');
ritzVector = ritzVector(:,order);

firstBlock = blockRanges{1};
greenCoupling = R0'*ritzVector(firstBlock,:);
normalization = real(trace(V'*V));
scalarWeights = real(sum(abs(greenCoupling).^2,1)).'/normalization;
assert(min(scalarWeights) >= -1e-13, ...
    'Block Lanczos produced a negative scalar spectral weight.');
scalarWeights(scalarWeights < 0 & ...
    scalarWeights > -100*eps(max(scalarWeights))) = 0;
scalarWeights = scalarWeights/sum(scalarWeights);

diagonalBlocks = cell(1,numel(blockRanges));
couplingBlocks = cell(1,max(numel(blockRanges)-1,0));
for block = 1:numel(blockRanges)
    diagonalBlocks{block} = ...
        blockJacobiMatrix(blockRanges{block},blockRanges{block});
    if block < numel(blockRanges)
        couplingBlocks{block} = blockJacobiMatrix( ...
            blockRanges{block+1},blockRanges{block});
    end
end

lanczos.basis = basis;
lanczos.projectedMatrix = blockJacobiMatrix;
lanczos.rawProjectedMatrix = rawProjectedMatrix;
lanczos.diagonalBlocks = diagonalBlocks;
lanczos.couplingBlocks = couplingBlocks;
lanczos.blockRanges = blockRanges;
lanczos.initialFactor = R0;
lanczos.ritzValues = ritzValue;
lanczos.ritzVectors = ritzVector;
lanczos.greenCoupling = greenCoupling;
lanczos.scalarWeights = scalarWeights;
lanczos.normalization = normalization;
lanczos.numberBlockSteps = numel(blocks);
lanczos.krylovDimension = size(basis,2);
lanczos.breakdown = breakdown;
lanczos.hermiticityError = hermiticityError;
lanczos.orthogonalityError = norm( ...
    basis'*basis-eye(size(basis,2)),'fro');
lanczos.blockTridiagonalDefect = norm( ...
    rawProjectedMatrix-blockJacobiMatrix,'fro') ...
    /max(norm(rawProjectedMatrix,'fro'),eps);
lanczos.options = options;

if options.verbose
    fprintf(['Block Lanczos: %d blocks, Krylov dimension %d, ' ...
        'orthogonality error %.3e, block defect %.3e.\n'], ...
        lanczos.numberBlockSteps,lanczos.krylovDimension, ...
        lanczos.orthogonalityError,lanczos.blockTridiagonalDefect);
end
end
