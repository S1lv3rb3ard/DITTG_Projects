function [relaxationFields,info] = getRelaxationFields(stack,dataFolder,options)
%GETRELAXATIONFIELDS Mechanically relaxed fields for a stack, with fallback.
%
% [relaxationFields,info] = getRelaxationFields(stack,dataFolder,options)
% looks in dataFolder for the Julia output of example2.jl whose twist
% angles match the stack,
%   triG_data_<t12>_<t23>_<N>.jld and triG_minimizer_<t12>_<t23>_<N>.jld,
% where t12 = theta_2 - theta_1 and t23 = theta_2 - theta_3 are printed with
% two decimals (the Julia rotation convention is clockwise). If several N
% are present, the largest is used unless options.N is set. The fields are
% loaded with loadJuliaTrilayerRelaxation.
%
% If no matching file exists, the behaviour is set by options.fallback:
%   'error' (default)  stop with instructions for running example2.jl
%   'toy'              use makeToyTrilayerRelaxation(stack,options.toyAmplitude)
%   'none'             zero relaxation
% info.source records which fields were returned.
if nargin < 3
    options = struct;
end
defaults = struct( ...
    'N',[], ...
    'fallback','error', ...
    'toyAmplitude',0.01*stack.aG, ...
    'loader',struct);
options = mergeOptions(options,defaults);

theta12 = stack.theta(2)-stack.theta(1);
theta23 = stack.theta(2)-stack.theta(3);
tag = sprintf('%.2f_%.2f',theta12,theta23);
pattern = fullfile(dataFolder,['triG_minimizer_',tag,'_*.jld']);
listing = dir(pattern);
sizes = zeros(1,numel(listing));
for index = 1:numel(listing)
    token = regexp(listing(index).name,'_(\d+)\.jld$','tokens','once');
    sizes(index) = str2double(token{1});
end
if ~isempty(options.N)
    listing = listing(sizes == options.N);
    sizes = sizes(sizes == options.N);
end

if ~isempty(listing)
    [N,best] = max(sizes);
    minimizerFile = fullfile(listing(best).folder,listing(best).name);
    dataFile = fullfile(listing(best).folder, ...
        sprintf('triG_data_%s_%d.jld',tag,N));
    [relaxationFields,info] = loadJuliaTrilayerRelaxation( ...
        minimizerFile,dataFile,stack,options.loader);
    info.source = 'julia';
    info.minimizerFile = minimizerFile;
    info.dataFile = dataFile;
    return
end

message = sprintf(['No Julia relaxation found for theta12 = %.2f, ' ...
    'theta23 = %.2f in %s. Run\n  julia example2.jl %.2f %.2f <N>\n' ...
    'and copy data/triG_*_%s_<N>.jld into that folder.'], ...
    theta12,theta23,dataFolder,theta12,theta23,tag);
switch lower(options.fallback)
    case 'error'
        error('getRelaxationFields:notFound','%s',message);
    case 'toy'
        warning('getRelaxationFields:toyFallback', ...
            '%s\nUsing the toy relaxation field instead.',message);
        relaxationFields = makeToyTrilayerRelaxation(stack,options.toyAmplitude);
        info.source = 'toy';
    case 'none'
        warning('getRelaxationFields:noRelaxation', ...
            '%s\nUsing zero relaxation instead.',message);
        zeroField = @(xk,xl) zeros(size(xk),'like',xk);
        relaxationFields = {zeroField,zeroField,zeroField};
        info.source = 'none';
    otherwise
        error('Unknown fallback "%s".',options.fallback);
end
end
