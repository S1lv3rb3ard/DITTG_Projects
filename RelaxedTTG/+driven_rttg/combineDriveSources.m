function sources = combineDriveSources(opticalSources,acousticSources)
%COMBINEDRIVESOURCES Concatenate possibly empty source structure arrays.
if isempty(opticalSources)
    sources = acousticSources(:);
elseif isempty(acousticSources)
    sources = opticalSources(:);
else
    sources = [opticalSources(:);acousticSources(:)];
end
end
