function [outputDir, cleanup] = experiment_setup(outputDir)
%EXPERIMENT_SETUP Resolve paths and restore caller settings on return.
rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if isempty(outputDir), outputDir = fullfile(rootDir,'results'); end
if ~isfolder(outputDir), mkdir(outputDir); end
state.rng = rng;
state.path = path;
state.text = get(groot,'defaultTextInterpreter');
state.ticks = get(groot,'defaultAxesTickLabelInterpreter');
state.legend = get(groot,'defaultLegendInterpreter');
cleanup = onCleanup(@() restore_settings(state));
addpath(fullfile(rootDir,'Tools'));
end

function restore_settings(state)
rng(state.rng);
path(state.path);
set(groot,'defaultTextInterpreter',state.text, ...
    'defaultAxesTickLabelInterpreter',state.ticks, ...
    'defaultLegendInterpreter',state.legend);
end
