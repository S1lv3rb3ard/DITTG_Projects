function projectRoot = setupRelaxedTTG()
%SETUPRELAXEDTTG Add the project root and non-package helpers to the path.
%
% Public observables are accessed through relaxed_dos.*, relaxed_idos.*,
% and driven_rttg.*. Shared numerical helpers use rttg_common.*.
projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(projectRoot));
end
