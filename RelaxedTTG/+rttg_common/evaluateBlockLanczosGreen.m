function [m,M] = evaluateBlockLanczosGreen(lanczos,z)
%EVALUATEBLOCKLANCZOSGREEN Evaluate V'*(H-z I)^(-1)*V after one recursion.
%
% m is the normalized trace compression. M(:,:,k) is the full block Green
% matrix at z(k). With the (H-z I)^(-1) convention, Im(m)>0 in the upper
% half-plane.
assert(isnumeric(z) && isvector(z) && ~isempty(z), ...
    'z must be a nonempty numeric vector.');
inputSize = size(z);
zRow = z(:).';
nodes = lanczos.ritzValues(:);
weights = lanczos.scalarWeights(:);

inverseDenominator = 1./(nodes-zRow);
m = weights.'*inverseDenominator;
m = reshape(m,inputSize);

if nargout > 1
    coupling = lanczos.greenCoupling;
    blockSize = size(coupling,1);
    M = complex(zeros(blockSize,blockSize,numel(zRow)));
    for index = 1:numel(zRow)
        M(:,:,index) = (coupling.*inverseDenominator(:,index).') ...
            *coupling';
    end
end
end
