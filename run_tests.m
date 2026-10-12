function results = run_tests
%RUN_TESTS Check group identities, numerical precision, and INS propagation.
rootDir = fileparts(mfilename('fullpath'));
results = runtests(fullfile(rootDir,'tests'));
assertSuccess(results);
end
