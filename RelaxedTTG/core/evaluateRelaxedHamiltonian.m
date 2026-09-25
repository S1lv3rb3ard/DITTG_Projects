function [H,parts,info] = evaluateRelaxedHamiltonian(cache,q)
%EVALUATERELAXEDHAMILTONIAN Assemble H(q) from a prepared relaxed cache.
assert(isequal(size(q),[2,1]),'q must be a 2-by-1 column vector.');
[Hintra,intralayerBlocks,intralayerInfo] = ...
    evaluateRelaxedIntralayerHamiltonian(cache.intralayer,q);
[Hinter,interlayerBlocks,interlayerInfo] = ...
    evaluateRelaxedInterlayerHamiltonian(cache.interlayer,q);
H = Hintra+Hinter;

relativeError = norm(H-H','fro')/max(norm(H,'fro'),eps);
assert(relativeError < cache.options.hermiticityTolerance, ...
    'The cached relaxed Hamiltonian is not Hermitian.');

if nargout > 1
    parts.intralayer = Hintra;
    parts.interlayer = Hinter;
    parts.intralayerBlocks = intralayerBlocks;
    parts.interlayerBlocks = interlayerBlocks;
end
if nargout > 2
    info.intralayer = intralayerInfo;
    info.interlayer = interlayerInfo;
    info.dimension = size(H,1);
    info.numNonzeros = nnz(H);
    info.relativeHermiticityError = relativeError;
end
end
