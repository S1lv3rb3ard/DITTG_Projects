function options = mergeOptions(options,defaults)
%MERGEOPTIONS Fill missing or empty fields from a defaults structure.
if nargin < 1 || isempty(options)
    options = struct;
end

names = fieldnames(defaults);
for k = 1:numel(names)
    name = names{k};
    if ~isfield(options,name) || isempty(options.(name))
        options.(name) = defaults.(name);
    elseif isstruct(options.(name)) && isstruct(defaults.(name))
        options.(name) = mergeOptions(options.(name),defaults.(name));
    end
end
end
