function DoF = getDoF(stack,W,L,output)
%GETDOF Generate the truncated reciprocal-space degrees of freedom.
if nargin < 4 || isempty(output)
    output = 'clean';
end
B = stack.B;

minSigma = min(cellfun(@(x) min(svd(x)),B));
mMax = ceil(L/minSigma);
[m1,m2] = ndgrid(-mMax:mMax);
m = [m1(:),m2(:)]';

minSigma = inf;
for layer = 1:3
    k = mod(layer,3)+1;
    l = mod(layer+1,3)+1;
    minSigma = min(minSigma,min(svd(B{k}-B{layer})));
    minSigma = min(minSigma,min(svd(B{l}-B{layer})));
end
nMax = ceil(W/minSigma)+1;
[n1,n2] = ndgrid(-nMax:nMax);
n = [n1(:),n2(:)]';

maxCount = 6*(2*mMax+1)^2*(2*nMax+1)^2;
DoF = zeros(maxCount,11);
shiftK = stack.K-stack.K(:,2);
count = 0;

for layer = 1:3
    k = mod(layer,3)+1;
    l = mod(layer+1,3)+1;
    normsK = vecnorm(B{k}*m,2,1);
    retainedM = normsK < L;
    C = B{k}-B{layer};
    D = B{l}-B{layer};
    xValues = m(:,retainedM);
    normsK = normsK(retainedM);
    centerMap = -(D\C);

    for column = 1:size(xValues,2)
        x = xValues(:,column);
        y = n+round(centerMap*x);
        normsL = vecnorm(B{l}*y,2,1);
        retainedL = normsL < L;
        y = y(:,retainedL);
        normsL = normsL(retainedL);
        if isempty(y)
            continue
        end

        residual = C*x+D*y;
        normsW = vecnorm(residual,2,1);
        retainedW = normsW < W;
        y = y(:,retainedW);
        residual = residual(:,retainedW);
        normsL = normsL(retainedW);
        normsW = normsW(retainedW);
        number = size(y,2);
        if number == 0
            continue
        end

        index = count+(1:number);
        count = count+number;
        DoF(index,1:2) = (residual+shiftK(:,layer)).';
        DoF(index,2*k+(1:2)) = repmat((B{k}*x).',number,1);
        DoF(index,2*l+(1:2)) = (B{l}*y).';
        DoF(index,9) = layer;
        DoF(index,10) = max(normsK(column),normsL).';
        DoF(index,11) = normsW.';
    end
end

DoF = sortrows(DoF(1:count,:),[9,1:8]);
if strcmpi(output,'clean')
    connectivity = getSparsity(DoF,'full');
    DoF = DoF(sum(connectivity,2) >= 2,:);
end
end
