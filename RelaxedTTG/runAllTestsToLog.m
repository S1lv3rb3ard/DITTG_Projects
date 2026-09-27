function results = runAllTestsToLog()
%RUNALLTESTSTOLOG Run every RelaxedTTG test and save the full report.
% Writes test_report.txt next to this file (command-window output of
% runtests, including every failure diagnostic, plus a summary table).
root = fileparts(mfilename('fullpath'));
addpath(genpath(root));
logFile = fullfile(root,'test_report.txt');
if exist(logFile,'file')
    delete(logFile);
end
diary(logFile);
cleanup = onCleanup(@() diary('off'));
fprintf('MATLAB %s, %s\n',version,datestr(now));
results = runtests(fullfile(root,'tests'),'IncludeSubfolders',true);
fprintf('\n===== SUMMARY =====\n');
disp(table(results));
fprintf('Passed %d, Failed %d, Incomplete %d, Total %d\n', ...
    nnz([results.Passed]),nnz([results.Failed]), ...
    nnz([results.Incomplete]),numel(results));
end
