function availableMemory = getAvailableMemory()
%GETAVAILABLEMEMORY Return approximately available system memory in bytes.
if ispc
    [~,systemView] = memory;
    availableMemory = systemView.PhysicalMemory.Available;
elseif ismac
    [status,text] = system('vm_stat');
    assert(status == 0,'Unable to query macOS memory.');
    pageSize = sscanf(text,'Mach Virtual Memory Statistics: (page size of %d bytes)');
    freePages = regexp(text,'Pages free:\s+(\d+)','tokens','once');
    inactivePages = regexp(text,'Pages inactive:\s+(\d+)','tokens','once');
    availableMemory = pageSize*(str2double(freePages{1}) ...
        +str2double(inactivePages{1}));
else
    [status,text] = system('awk ''/MemAvailable/ {print $2}'' /proc/meminfo');
    assert(status == 0,'Unable to query Linux memory.');
    availableMemory = 1024*str2double(strtrim(text));
end
end
