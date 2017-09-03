function run_test( varargin )
%RUN_TEST Verify all the tests.

addpath(pwd); % add the pwd path to the environment

test_ndg_lib = 1;

if test_ndg_lib
    filepath{1} = 'testing/StdRegions/Triangle';
    filepath{2} = 'testing/StdRegions/Quad';
    test_dir(filepath);
end

end

function test_dir(filepath)
% test all the tests in the specific folder
for f = 1:numel(filepath) % loop over all the path
    file = dir(filepath{f});
    for i = 1:numel(file) % loop over all the files
        [~, ~, ext] = fileparts(file(i).name); 
        if strcmp(ext, '.m') % test the '.m' files
            results = runtests(fullfile(filepath{f}, file(i).name));
            table(results)
        end
    end
end
end% func