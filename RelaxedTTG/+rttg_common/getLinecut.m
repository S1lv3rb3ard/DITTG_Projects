function linecut = getLinecut(K,Nq)
%GETLINECUT Preserve the paper's K1-K2-Gamma12-Gamma23-K2-K3 path.
midpoint = (K(:,1:2)+K(:,2:3))/2;
edge = K(:,2:3)-K(:,1:2);
normal = rotationMatrixDegrees(-90)*edge;
normal = normal./vecnorm(normal,2,1);
center = midpoint+normal.*vecnorm(edge,2,1);
vertices = [K(:,1:2),center,K(:,2:3)];

linecut.q = vertices(:,1);
linecut.tickLocs = zeros(1,size(vertices,2));
linecut.kPath = 0;
for segment = 1:size(vertices,2)-1
    startPoint = vertices(:,segment);
    endPoint = vertices(:,segment+1);
    distance = norm(endPoint-startPoint);
    parameter = linspace(0,1,Nq+1);
    points = startPoint.*(1-parameter)+endPoint.*parameter;
    path = linecut.tickLocs(segment)+linspace(0,distance,Nq+1);
    linecut.q = [linecut.q,points(:,2:end)]; %#ok<AGROW>
    linecut.kPath = [linecut.kPath,path(2:end)]; %#ok<AGROW>
    linecut.tickLocs(segment+1) = linecut.tickLocs(segment)+distance;
end
linecut.labels = {'$K_{1}$','$K_{2}$','$\Gamma_{12}$', ...
    '$\Gamma_{23}$','$K_{2}$','$K_{3}$'};
end
