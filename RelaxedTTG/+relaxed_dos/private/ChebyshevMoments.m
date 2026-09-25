function moments = ChebyshevMoments(H,Y,P,mode,precision)
%CHEBYSHEVMOMENTS Compute projected KPM moments for a batch of q-points.
if nargin < 5 || isempty(precision)
    precision = 'double';
end
number = numel(H);
vectorLength = cellfun(@(matrix) size(matrix,1),H);
sumByPoint = arrayfun(@(n) ones(1,n),vectorLength, ...
    'UniformOutput',false);
sumByPoint = sparse(blkdiag(sumByPoint{:}));
H = sparse(blkdiag(H{:}));

selectionSize = cellfun(@(matrix) size(matrix,2),Y);
maximumSelectionSize = max(selectionSize);
for index = 1:number
    if selectionSize(index) < maximumSelectionSize
        Y{index} = [Y{index},zeros(size(Y{index},1), ...
            maximumSelectionSize-selectionSize(index))];
    end
end
Y0 = cat(1,Y{:});

if strcmpi(precision,'single')
    H = single(H);
    sumByPoint = single(sumByPoint);
    Y0 = single(Y0);
end
moments = zeros(number,P+1);
if strcmpi(mode,'gpu')
    H = gpuArray(H);
    sumByPoint = gpuArray(sumByPoint);
    Y0 = gpuArray(Y0);
    moments = gpuArray(moments);
end
Y0 = full(Y0);
innerProduct = @(A,B) sumByPoint*sum(real(conj(A).*B),2);

previous = Y0;
current = H*Y0;
mu0 = innerProduct(previous,previous);
mu1 = innerProduct(current,previous);
moments(:,1) = mu0;
if P >= 1
    moments(:,2) = mu1;
end

H = sparse(2.*H);
for order = 1:floor(P/2)
    moments(:,2*order+1) = 2*innerProduct(current,current)-mu0;
    next = H*current-previous;
    if 2*order+2 <= P+1
        moments(:,2*order+2) = ...
            2*innerProduct(next,current)-mu1;
    end
    previous = current;
    current = next;
end
end
