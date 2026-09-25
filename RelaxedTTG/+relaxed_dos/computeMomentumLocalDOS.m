function Job = computeMomentumLocalDOS(Job,E,mode)
%COMPUTEMOMENTUMLOCALDOS Paper-defined relaxed non-driven momentum LDoS.
%
% This is the public module entry point for equations (3.19), (3.20), and
% (4.51). It averages the six (G=0,j,alpha) diagonal spectral entries.
Job = getLDoS(Job,E,mode);
end
