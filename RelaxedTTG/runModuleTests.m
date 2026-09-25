function results = runModuleTests(moduleName)
%RUNMODULETESTS Run tests for exactly one public observable module.
arguments
    moduleName (1,1) string {mustBeMember(moduleName, ...
        ["relaxed_dos","relaxed_idos","driven_rttg","common","all"])}
end

projectRoot = fileparts(mfilename('fullpath'));
switch moduleName
    case "all"
        target = fullfile(projectRoot,'tests');
    otherwise
        target = fullfile(projectRoot,'tests',moduleName);
end
results = runtests(target,'IncludeSubfolders',true);
end
