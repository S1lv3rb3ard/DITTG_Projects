function selector = getDrivenCentralSelector(driveCache,spatialSelector)
%GETDRIVENCENTRALSELECTOR Embed a spatial selector in temporal sector zero.
assert(size(spatialSelector,1) == driveCache.baseDimension, ...
    'spatialSelector has the wrong number of rows.');
selector = sparse(driveCache.extendedDimension,size(spatialSelector,2));
selector(driveCache.centralSpatialRange,:) = spatialSelector;
end
