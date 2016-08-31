function configure

%% Polylib
localpath = '+Polylib';
installpath = localpath;
src = {'JacobiP.c', ...
    'zwglj.c', ...
    'Dglj.c', ...
    'jacobfd.c', ...
    'GradJacobiP.c'};
libsrc = {'polylib.c'};
fprintf('==============Polylib=================\n')
install(localpath, installpath, src, libsrc);
fprintf('==============Polylib=================\n\n')

%% Limiter1d
localpath = '+Utilities/+Limiter/Mex';
installpath = '+Utilities/+Limiter/+Limiter1D';
src = {'Minmod1d_Mex.c'};
libsrc = {'MatrixUtilities.c', 'BasicFunction.c'};
fprintf('==============Limiter2d=================\n')
install(localpath, installpath, src, libsrc);
fprintf('==============Limiter2d=================\n\n')


%% Limiter2d
localpath = '+Utilities/+Limiter/Mex';
installpath = '+Utilities/+Limiter/+Limiter2D';
src = {'BJ2d_Mex.c', ...
    'BJLoc2d_Mex.c',...
    'HWENO2d_Mex.c',...
    'KXRCF_detector2d_Mex.c',...
    'TVB2d_Mex.c',...
    'VB2d_Mex.c',...
    'TVB_detector2d_Mex.c'};
libsrc = {'MatrixUtilities.c', 'BasicFunction.c'};
fprintf('==============Limiter2d=================\n')
install(localpath, installpath, src, libsrc);
fprintf('==============Limiter2d=================\n\n')

%% SWE2d

end% func

%% install
% Compile the source file and put into spicific directory.
function install(localpath, installpath, src, libsrc)
% 
pwdPath  = pwd;
fullPath = fullfile(pwd, localpath);
dirPath  = fullfile(pwd, installpath);

srcfile  = fullfile(fullPath, src);
libfile  = fullfile(fullPath, libsrc);

cd(fullPath);

for i = 1:numel(srcfile)
    fprintf('installing %s to %s...\n', src{i}, installpath);
    file = [srcfile(i), libfile{:}];
    mex('-O', file{:}, '-outdir', dirPath);
end% for

cd(pwdPath);
end