localPath = pwd;
polylibPath = fullfile(pwd, '+Polylib');
cd(polylibPath);
fprintf('=================================\n')
fprintf('Start installing Polylib\n')
fprintf('Installing of JacobiP.c\n')
mex('-O', [polylibPath,'/JacobiP.c'], [polylibPath, '/polylib.c'])
fprintf('Installing of zwglj.c\n')
mex('-O', [polylibPath, '/zwglj.c'], [polylibPath, '/polylib.c'])
fprintf('Installing of Dglj.c\n')
mex('-O', [polylibPath, '/Dglj.c'], [polylibPath, '/polylib.c'])
fprintf('Installing of jacobfd.c\n')
mex('-O', [polylibPath, '/jacobfd.c'], [polylibPath, '/polylib.c'])
fprintf('Installing of GradJacobiP.c\n')
mex('-O', [polylibPath, '/GradJacobiP.c'], [polylibPath, '/polylib.c'])
fprintf('=================================\n')
fprintf('Finish installation Polylib\n')
cd(localPath);