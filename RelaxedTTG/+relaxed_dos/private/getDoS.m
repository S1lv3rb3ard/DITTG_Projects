function Job = getDoS(Job,E,mode)
%GETDOS Reconstruct total DoS from KPM moments and quadrature weights.
g = Jackson(Job.P);
thetaValues = sort(unique(Job.list.t)).';
for thetaIndex = thetaValues
    thetaRows = Job.list.t == thetaIndex;
    scale = Job.scale(thetaIndex);
    scaledE = scale*E;
    T = Chebyshev(Job.P,scaledE,mode);
    T(1,:) = 0.5*T(1,:);
    ChebyshevWeight = (2*scale)./(pi*sqrt(1-scaledE.^2));
    for selectorIndex = 1:3
        rows = thetaRows & Job.list.X == selectorIndex;
        N = Job.list.N(rows);
        weight = Job.weights{thetaIndex}(N,selectorIndex);
        Job.list.DoS(rows,:) = weight ...
            .*((g.*Job.list.mu(rows,:))*T).*ChebyshevWeight;
    end
end
Job.list.DoS = gather(Job.list.DoS);
Job.list = groupsummary(Job.list,["t","N","model"],"sum","DoS");
Job.list.GroupCount = [];
Job.list = renamevars(Job.list,'sum_DoS','DoS');
end
