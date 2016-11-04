function runTest( varargin )
%RUNTEST 测试脚本
% 

% 为了使测试函数顺利执行，必须将 NDG-FEM 库路径添加到环境变量内
addpath(pwd);

if nargin == 0
    test_StdRegions = true;
else
    test_StdRegions = false;
    
    for i = nargin
        switch varargin{i}
            case 'StdRegions'
                test_StdRegions = true;
            
            otherwise
                error(['No test for %s, please choose one of\n',...
                    '  StdRegions\n'], varargin{i})
        end% switch
    end% for
end% if

if test_StdRegions
    filepath{1} = 'testing/StdRegions/Triangle';
    testDir(filepath);
end

end

function testDir(filepath)
% TESTDIR 测试给定路径下所有测试脚本
for f = 1:numel(filepath) % 遍历所有文件路径
    file = dir(filepath{f});
    for i = 1:numel(file) % 遍历路径下所有文件
        [~, ~, ext] = fileparts(file(i).name); 
        if strcmp(ext, '.m') % 文件后缀为.m则执行测试
            results = runtests(fullfile(filepath, file(i).name));
            table(results)
        end
    end
end
end% func