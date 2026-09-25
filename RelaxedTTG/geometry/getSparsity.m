function S = getSparsity(DoF,output)
%GETSPARSITY Unrelaxed nearest-layer connectivity of the reciprocal basis.
intraS = logical(speye(size(DoF,1)));
I = cell(1,2);
J = cell(1,2);

for pair = 1:2
    rows = find(DoF(:,9) == pair);
    columns = find(DoF(:,9) == pair+1);
    spectatorColumns = 2*(mod(pair+1,3)+1)+(1:2);
    X = DoF(rows,spectatorColumns);
    Y = DoF(columns,spectatorColumns);
    [uniqueY,~,groupY] = unique(Y,'rows');
    groupY = accumarray(groupY,(1:numel(columns))',[],@(x){x});

    [matched,location] = ismember(X,uniqueY,'rows');
    matchedCells = groupY(location(matched));
    I{pair} = repelem(rows(matched),cellfun(@numel,matchedCells));
    J{pair} = columns(vertcat(matchedCells{:}));
end

finalI = [I{1};I{2};J{1};J{2}];
finalJ = [J{1};J{2};I{1};I{2}];
interS = sparse(finalI,finalJ,true,size(DoF,1),size(DoF,1));
S = intraS | interS;

if strcmpi(output,'block')
    numDoF = getNumDoF(DoF);
    S = mat2cell(S,numDoF,numDoF);
end
end
