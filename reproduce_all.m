function results = reproduce_all(outputDir)
%REPRODUCE_ALL Generate the numerical examples and figures from the tutorial.
%   REPRODUCE_ALL writes to results/ beside this file. An optional directory
%   argument selects another destination. Existing generated files there
%   are replaced. All published Monte Carlo realizations are run.

rootDir = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(outputDir)
    outputDir = fullfile(rootDir,'results');
end
oldPath = path;
cleanup = onCleanup(@() path(oldPath));
addpath(fullfile(rootDir,'Examples'));
if ~isfolder(outputDir), mkdir(outputDir); end
started = tic;
fprintf('Reproducing tutorial figures in %s\n', outputDir);
results.discretization = discretization(outputDir);
results.loglinear = loglinear_error(outputDir);
results.distributions = error_distributions(outputDir);
results.uncertainty = uncertainty_propagation(outputDir);
results.filters = compare_filters(outputDir);
results.elapsed_seconds = toc(started);
fprintf('All examples completed in %.1f seconds.\n',results.elapsed_seconds);
end
