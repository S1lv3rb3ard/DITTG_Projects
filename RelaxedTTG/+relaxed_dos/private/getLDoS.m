function Job = getLDoS(Job,E,mode)
%GETLDOS Reconstruct the paper's averaged momentum LDoS from KPM moments.
%
% With Job.X selecting all six (0,j,alpha) orbitals, the factor 1/6 is
% exactly 1/|A| in equations (3.20) and (4.51).
g = Jackson(Job.P);
thetaValues = sort(unique(Job.list.t)).';
for thetaIndex = thetaValues
    rows = Job.list.t == thetaIndex;
    scale = Job.scale(thetaIndex);
    scaledE = scale*E;
    T = Chebyshev(Job.P,scaledE,mode);
    T(1,:) = 0.5*T(1,:);
    ChebyshevWeight = (2*scale)./(pi*sqrt(1-scaledE.^2));
    Job.list.DoS(rows,:) = (1/6) ...
        *((g.*Job.list.mu(rows,:))*T).*ChebyshevWeight;
end
Job.list.DoS = toHost(Job.list.DoS,strcmpi(mode,'gpu'));
end
