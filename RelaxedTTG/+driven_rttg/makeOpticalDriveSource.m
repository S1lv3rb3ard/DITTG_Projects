function source = makeOpticalDriveSource( ...
    amplitude,waveVector,frequency,phase,cutoff,envelope,couplingOperator)
%MAKEOPTICALDRIVESOURCE Define one pulsed traveling optical source.
%
% amplitude is an energy, waveVector is 2-by-1, and frequency is an angular
% frequency. The Hamiltonian energy shift is options.hbar*frequency. Pulse
% times therefore use the reciprocal unit of frequency.
if nargin < 7
    couplingOperator = [];
end
source.kind = 'optical';
source.amplitude = amplitude;
source.waveVector = waveVector(:);
source.frequency = frequency;
source.phase = phase;
source.cutoff = cutoff;
source.envelope = envelope;
source.couplingOperator = couplingOperator;
end
