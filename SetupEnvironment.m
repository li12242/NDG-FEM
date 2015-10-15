addpath(pwd)
addpath([pwd, '/Codes1D'], [pwd, '/Codes1D/testing'])
addpath([pwd, '/Codes2D'], [pwd, '/Codes2D/testing'])
addpath([pwd, '/Codes3D'], [pwd, '/Codes3D/testing'])
addpath([pwd, '/ServiceRoutines/'], [pwd, '/Grid'])
addpath([pwd, '/src'])
addpath([pwd, '/src/testing'])


files = dir('Grid');
for i = 4: numel(files)
    addpath([pwd, '/Grid/', files(i).name])
end