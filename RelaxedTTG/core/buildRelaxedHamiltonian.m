function [H,parts,info] = buildRelaxedHamiltonian( ...
    stack,DoF,q,shellsReference,relaxationFields,options)
%BUILDRELAXEDHAMILTONIAN Assemble the full relaxed trilayer Hamiltonian.
%
% H = H_intralayer + H_interlayer in the common 2*size(DoF,1) basis.
% options.intralayer and options.interlayer are forwarded to the two
% component builders.
if nargin < 6
    options = struct;
end
defaults = struct( ...
    'intralayer',struct, ...
    'interlayer',struct, ...
    'hermiticityTolerance',1e-12, ...
    'verbose',true);
options = mergeOptions(options,defaults);

[Hintra,intralayerBlocks,intralayerInfo] = ...
    buildRelaxedIntralayerHamiltonian( ...
        stack,DoF,q,shellsReference,relaxationFields,options.intralayer);
[Hinter,interlayerBlocks,interlayerInfo] = ...
    buildRelaxedInterlayerHamiltonian( ...
        stack,DoF,q,relaxationFields,options.interlayer);

H = Hintra+Hinter;
relativeError = norm(H-H','fro')/max(norm(H,'fro'),eps);
assert(relativeError < options.hermiticityTolerance, ...
    'The full relaxed Hamiltonian is not Hermitian.');

parts.intralayer = Hintra;
parts.interlayer = Hinter;
parts.intralayerBlocks = intralayerBlocks;
parts.interlayerBlocks = interlayerBlocks;

info.intralayer = intralayerInfo;
info.interlayer = interlayerInfo;
info.dimension = size(H,1);
info.numNonzeros = nnz(H);
info.relativeHermiticityError = relativeError;
info.options = options;

if options.verbose
    fprintf(['Full relaxed Hamiltonian: %d x %d, nnz=%d, ' ...
        'relative Hermiticity error %.3e.\n'], ...
        size(H,1),size(H,2),nnz(H),relativeError);
end
end
