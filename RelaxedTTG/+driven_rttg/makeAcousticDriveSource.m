function source = makeAcousticDriveSource( ...
    amplitude,waveVector,frequency,phase,cutoff,couplingOperator,envelope)
%MAKEACOUSTICDRIVESOURCE Define one coherent traveling acoustic source.
%
% The default coupling is the scalar deformation potential of the note.
% couplingOperator may instead contain a layer-orbital deformation matrix.
if nargin < 6
    couplingOperator = [];
end
if nargin < 7 || isempty(envelope)
    envelope = struct('type','constant');
end
source.kind = 'acoustic';
source.amplitude = amplitude;
source.waveVector = waveVector(:);
source.frequency = frequency;
source.phase = phase;
source.cutoff = cutoff;
source.envelope = envelope;
source.couplingOperator = couplingOperator;
end
