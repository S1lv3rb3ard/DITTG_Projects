function [Hintra,blocks,info] = evaluateRelaxedIntralayerHamiltonian(cache,q)
%EVALUATERELAXEDINTRALAYERHAMILTONIAN Evaluate cached intralayer terms.
assert(isequal(size(q),[2,1]),'q must be a 2-by-1 column vector.');

Hintra = sparse(cache.dimension,cache.dimension);
blocks = cell(1,3);
info = cell(1,3);
for layer = 1:3
    data = cache.layer{layer};
    phaseA = exp(-1i*q.'*data.RA).';
    phaseAB = exp(-1i*q.'*data.RBminusTau).';
    phaseBA = exp(-1i*q.'*data.RBplusTau).';

    values = [data.coeffAA*phaseA, ...
              data.coeffAB*phaseAB, ...
              data.coeffBA*phaseBA, ...
              data.coeffBB*phaseA];
    values = data.couplingChi.*values;

    diagonal = data.diagonal;
    if any(diagonal)
        aa = values(diagonal,1);
        ab = values(diagonal,2);
        ba = values(diagonal,3);
        bb = values(diagonal,4);
        values(diagonal,:) = [0.5*real(aa), ...
            0.25*(ab+conj(ba)),0.25*(ba+conj(ab)),0.5*real(bb)];
    end

    number = numel(data.source);
    rowIndex = zeros(4*number,1);
    colIndex = zeros(4*number,1);
    flatValues = complex(zeros(4*number,1));
    rowBase = 2*(data.source-1);
    colBase = 2*(data.target-1);
    rowIndex(1:4:end) = rowBase+1;
    rowIndex(2:4:end) = rowBase+1;
    rowIndex(3:4:end) = rowBase+2;
    rowIndex(4:4:end) = rowBase+2;
    colIndex(1:4:end) = colBase+1;
    colIndex(2:4:end) = colBase+2;
    colIndex(3:4:end) = colBase+1;
    colIndex(4:4:end) = colBase+2;
    flatValues(1:4:end) = values(:,1);
    flatValues(2:4:end) = values(:,2);
    flatValues(3:4:end) = values(:,3);
    flatValues(4:4:end) = values(:,4);

    upper = sparse(rowIndex,colIndex,flatValues, ...
        2*data.numStates,2*data.numStates);
    block = upper+upper';
    blocks{layer} = block;
    orbitalIndex = reshape( ...
        [2*data.dofIndex-1,2*data.dofIndex].',[],1);
    Hintra(orbitalIndex,orbitalIndex) = block;

    info{layer} = struct( ...
        'layer',layer, ...
        'dofIndex',data.dofIndex, ...
        'numStates',data.numStates, ...
        'numRetainedUnorderedPairs',number, ...
        'numNonzeros',nnz(block));
end
end
