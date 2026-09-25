function Job = BatchWeights(Job,mode)
%BATCHWEIGHTS Compute KPM moments for every row of Job.list.
[BatchIdx,BatchCount,numGPUs,~] = ...
    getBatches(Job.theta,Job.W,Job.L,Job.list,mode);

BatchList = cell(BatchCount,1);
for batch = 1:BatchCount
    BatchList{batch} = Job.list(BatchIdx{batch},:);
end

fprintf('Starting DoS computation: %d batches total...\n',BatchCount)
for batch = 1:BatchCount
    if strcmpi(mode,'gpu') && numGPUs > 0
        gpuDevice(1);
    end
    rows = BatchList{batch};
    nJobs = height(rows);
    cellH = cell(1,nJobs);
    selection = cell(1,nJobs);
    for job = 1:nJobs
        thetaIndex = rows.t(job);
        selectorIndex = rows.X(job);
        cellH{job} = Job.scale(thetaIndex) ...
            .*Job.H{thetaIndex}(rows.Q(job,:).');
        selection{job} = Job.X{thetaIndex}{selectorIndex};
    end
    rows.mu = toHost(ChebyshevMoments( ...
        cellH,selection,Job.P,mode), ...
        strcmpi(mode,'gpu'));
    BatchList{batch} = rows;
    fprintf('Computed: batch %d of %d\n',batch,BatchCount);
end

JobList = cell2mat(BatchList);
if ismember('cost',JobList.Properties.VariableNames)
    JobList.cost = [];
end
Job.list = JobList;
end
