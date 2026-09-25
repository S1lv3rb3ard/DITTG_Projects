function [BatchIdx,BatchCount,numGPUs,numCPUWorkers] = ...
    getBatches(theta,W,L,JobList,mode)
%GETBATCHES Pack KPM jobs according to the supplied memory-cost estimate.
memoryFraction = 0.5;
numberJobs = height(JobList);
numGPUs = 0;
numCPUWorkers = 0;

if strcmpi(mode,'gpu')
    numGPUs = gpuDeviceCount('available');
    assert(numGPUs > 0,'GPU mode requested but no GPUs are available.');
    availableMemory = inf;
    for device = 1:numGPUs
        availableMemory = min(availableMemory,gpuDevice(device).AvailableMemory);
    end
    availableMemory = memoryFraction*availableMemory;
elseif strcmpi(mode,'cpu')
    availableMemory = memoryFraction*getAvailableMemory();
else
    error('Computation mode must be gpu or cpu.');
end

badJob = find(JobList.cost > availableMemory,1);
if ~isempty(badJob)
    thetaIndex = JobList.t(badJob);
    error(['Job %d (W=%.2f, L=%g, theta index %d) requires %.2f GB, ' ...
        'exceeding the %.2f GB batch limit.'],badJob,W,L(thetaIndex), ...
        thetaIndex,JobList.cost(badJob)/1e9,availableMemory/1e9);
end

BatchIdx = cell(numberJobs,1);
BatchCount = 0;
first = 1;
currentMemory = 0;
for job = 1:numberJobs
    if currentMemory+JobList.cost(job) > availableMemory
        BatchCount = BatchCount+1;
        BatchIdx{BatchCount} = first:job-1;
        first = job;
        currentMemory = JobList.cost(job);
    else
        currentMemory = currentMemory+JobList.cost(job);
    end
end
if first <= numberJobs
    BatchCount = BatchCount+1;
    BatchIdx{BatchCount} = first:numberJobs;
end
BatchIdx = BatchIdx(1:BatchCount);
end
