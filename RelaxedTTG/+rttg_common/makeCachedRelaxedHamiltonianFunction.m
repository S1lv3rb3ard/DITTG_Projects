function Hfun = makeCachedRelaxedHamiltonianFunction(cache,maxMatrices)
%MAKECACHEDRELAXEDHAMILTONIANFUNCTION Return a memoized H(q) function.
%
% Repeated Job.list rows for different observables commonly request the same
% q.  This small exact-q LRU cache avoids assembling that matrix repeatedly.
if nargin < 2 || isempty(maxMatrices)
    maxMatrices = 1;
end
assert(maxMatrices >= 0 && maxMatrices == fix(maxMatrices), ...
    'maxMatrices must be a nonnegative integer.');

keys = cell(1,0);
matrices = cell(1,0);
Hfun = @evaluate;

    function H = evaluate(q)
        assert(isequal(size(q),[2,1]),'q must be a 2-by-1 column vector.');
        key = sprintf('%.17g,%.17g',q(1),q(2));
        location = find(strcmp(keys,key),1);
        if ~isempty(location)
            H = matrices{location};
            if location > 1
                keys = [{keys{location}},keys(1:location-1), ...
                    keys(location+1:end)];
                matrices = [{matrices{location}},matrices(1:location-1), ...
                    matrices(location+1:end)];
            end
            return
        end

        H = evaluateRelaxedHamiltonian(cache,q);
        if maxMatrices > 0
            keys = [{key},keys];
            matrices = [{H},matrices];
            keys = keys(1:min(maxMatrices,numel(keys)));
            matrices = matrices(1:min(maxMatrices,numel(matrices)));
        end
    end
end
