function Job = computeTotalMomentumDOS(Job,E,mode)
%COMPUTETOTALMOMENTUMDOS Paper-defined total relaxed non-driven DoS.
%
% Job must already contain the reciprocal-space weights, selectors, and
% Chebyshev moments produced by the existing paper DoS workflow.
Job = getDoS(Job,E,mode);
end
