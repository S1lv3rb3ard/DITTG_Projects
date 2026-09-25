function Job = computeMoments(Job,mode)
%COMPUTEMOMENTS Compute Jackson-Chebyshev moments for a prepared DoS Job.
%
% This public entry point owns the batching step shared by the momentum
% local and total DoS reconstructions in this module.
Job = BatchWeights(Job,mode);
end
