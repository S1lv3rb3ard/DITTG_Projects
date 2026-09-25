function value = evaluateDriveEnvelope(source,t)
%EVALUATEDRIVEENVELOPE Evaluate a real scalar source envelope.
assert(isfield(source,'envelope') && isfield(source.envelope,'type'), ...
    'Every drive source must define envelope.type.');
envelope = source.envelope;
type = lower(string(envelope.type));
tShape = size(t);
t = t(:).';

switch type
    case "constant"
        value = ones(size(t));

    case "gaussiantrain"
        requireFields(envelope,{'period','fwhm','center'});
        assert(envelope.period > 0 && envelope.fwhm > 0 && ...
            envelope.fwhm < envelope.period, ...
            'A Gaussian train requires 0 < fwhm < period.');
        distance = wrappedDistance(t-envelope.center,envelope.period);
        value = exp(-4*log(2)*(distance/envelope.fwhm).^2);

    case "sin2train"
        requireFields(envelope,{'period','duration','center'});
        assert(envelope.period > 0 && envelope.duration > 0 && ...
            envelope.duration < envelope.period, ...
            'A sin2 train requires 0 < duration < period.');
        distance = wrappedDistance(t-envelope.center,envelope.period);
        inside = abs(distance) <= envelope.duration/2;
        value = zeros(size(t));
        value(inside) = cos(pi*distance(inside)/envelope.duration).^2;

    case "gaussian"
        requireFields(envelope,{'fwhm','centers'});
        assert(envelope.fwhm > 0 && ~isempty(envelope.centers), ...
            'A Gaussian pulse requires positive fwhm and pulse centers.');
        value = zeros(size(t));
        for center = envelope.centers(:).'
            value = value+exp(-4*log(2)*((t-center)/envelope.fwhm).^2);
        end
        value = min(value,1);

    case "custom"
        assert(isfield(envelope,'function') && ...
            isa(envelope.function,'function_handle'), ...
            'A custom envelope requires envelope.function.');
        value = envelope.function(t);

    otherwise
        error('Unknown envelope type "%s".',envelope.type);
end

assert(numel(value) == numel(t) && isreal(value) && ...
    all(isfinite(value),'all') && all(value >= 0,'all') && ...
    all(value <= 1+1e-12,'all'), ...
    'Drive envelopes must be real, finite, and between zero and one.');
value = reshape(value,tShape);
end

function distance = wrappedDistance(offset,period)
distance = mod(offset+period/2,period)-period/2;
end

function requireFields(structure,names)
for index = 1:numel(names)
    assert(isfield(structure,names{index}), ...
        'Envelope is missing field "%s".',names{index});
end
end
