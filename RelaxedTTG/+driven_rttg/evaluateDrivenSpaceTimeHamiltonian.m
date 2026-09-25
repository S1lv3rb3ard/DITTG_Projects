function [K,parts,info] = evaluateDrivenSpaceTimeHamiltonian(cache,t)
%EVALUATEDRIVENSPACETIMEHAMILTONIAN Assemble K_q(t) at envelope time t.
assert(isscalar(t) && isfinite(t),'t must be a finite scalar.');
K = cache.staticMatrix;
parts.static = cache.staticMatrix;
parts.optical = sparse(cache.extendedDimension,cache.extendedDimension);
parts.acoustic = sparse(cache.extendedDimension,cache.extendedDimension);
parts.source = cell(numel(cache.sources),1);
envelopeValue = zeros(numel(cache.sources),1);

for source = 1:numel(cache.sources)
    envelopeValue(source) = driven_rttg.evaluateDriveEnvelope( ...
        cache.sources(source),t);
    contribution = envelopeValue(source)*cache.sourceMatrices{source};
    parts.source{source} = contribution;
    K = K+contribution;
    if source <= numel(cache.opticalSources)
        parts.optical = parts.optical+contribution;
    else
        parts.acoustic = parts.acoustic+contribution;
    end
end

relativeHermiticityError = norm(K-K','fro')/max(norm(K,'fro'),eps);
assert(relativeHermiticityError <= cache.options.hermiticityTolerance, ...
    'The driven extended Hamiltonian is not Hermitian.');
info.time = t;
info.envelopeValue = envelopeValue;
info.relativeHermiticityError = relativeHermiticityError;
info.dimension = size(K,1);
info.numNonzeros = nnz(K);
end
