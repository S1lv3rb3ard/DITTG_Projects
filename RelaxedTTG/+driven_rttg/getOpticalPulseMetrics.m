function metrics = getOpticalPulseMetrics(source)
%GETOPTICALPULSEMETRICS Return intensity-equivalent pulse duty information.
%
% The equivalent duty cycle is the period average of envelope(t)^2 because
% optical intensity scales as the squared field envelope.
assert(strcmpi(source.kind,'optical'), ...
    'Pulse metrics apply only to optical sources.');
envelope = source.envelope;
type = lower(string(envelope.type));

switch type
    case "gaussiantrain"
        metrics.equivalentIntensityDutyCycle = ...
            sqrt(pi/(8*log(2)))*envelope.fwhm/envelope.period;
    case "sin2train"
        metrics.equivalentIntensityDutyCycle = ...
            3*envelope.duration/(8*envelope.period);
    case "gaussian"
        metrics.equivalentIntensityDutyCycle = NaN;
    case "custom"
        assert(isfield(envelope,'equivalentIntensityDutyCycle'), ...
            ['A custom optical envelope must provide ' ...
             'equivalentIntensityDutyCycle.']);
        metrics.equivalentIntensityDutyCycle = ...
            envelope.equivalentIntensityDutyCycle;
    otherwise
        metrics.equivalentIntensityDutyCycle = 1;
end

metrics.peakIntensity = NaN;
metrics.estimatedAverageIntensity = NaN;
metrics.maximumAverageIntensity = NaN;
if isfield(envelope,'peakIntensity')
    metrics.peakIntensity = envelope.peakIntensity;
    metrics.estimatedAverageIntensity = envelope.peakIntensity ...
        *metrics.equivalentIntensityDutyCycle;
end
if isfield(envelope,'maximumAverageIntensity')
    metrics.maximumAverageIntensity = envelope.maximumAverageIntensity;
end
end
