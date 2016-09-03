function configure(varargin)

if nargin == 0
    conf_Polylib = true;
    conf_Limit1d = true;
    conf_Limit2d = true;
    conf_SWE2d   = true;
else
    conf_Polylib = false;
    conf_Limit1d = false;
    conf_Limit2d = false;
    conf_SWE2d   = false;
    for i = nargin
        switch varargin{i}
            case 'Polylib'
                conf_Polylib = true;
            case 'Limiter1d'
                conf_Limit1d = true;
            case 'Limiter2d'
                conf_Limit2d = true;
            case 'SWE2d'
                conf_SWE2d = true;
            otherwise
                error(['No options for %s, please choose one of\n'...
                    ,'  Polylib\n  Limiter\n  SWE2d\n'], varargin{i})
        end% switch
    end% for
end% if

%% Polylib
if conf_Polylib
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
end
%% Limiter1d
if conf_Limit1d
    localpath = '+Utilities/+Limiter/Mex';
    installpath = '+Utilities/+Limiter/+Limiter1D';
    src = {'Minmod1d_Mex.c'};
    libsrc = {'MatrixUtilities.c', 'BasicFunction.c'};
    fprintf('==============Limiter2d=================\n')
    install(localpath, installpath, src, libsrc);
    fprintf('==============Limiter2d=================\n\n')
end

%% Limiter2d
if conf_Limit2d
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
end
%% SWE2d
if conf_SWE2d
end

end% func

%% install function
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