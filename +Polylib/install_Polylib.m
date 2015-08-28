localPath = pwd;
polylibPath = fullfile(pwd, '+Polylib');
cd(polylibPath);
mex('-O', [polylibPath,'/JacobiP.c'], [polylibPath, '/polylib.c'])
mex('-O', [polylibPath, '/zwglj.c'], [polylibPath, '/polylib.c'])
mex('-O', [polylibPath, '/Dglj.c'], [polylibPath, '/polylib.c'])
mex('-O', [polylibPath, '/jacobfd.c'], [polylibPath, '/polylib.c'])
mex('-O', [polylibPath, '/GradJacobiP.c'], [polylibPath, '/polylib.c'])
cd(localPath);