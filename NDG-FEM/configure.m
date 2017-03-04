function configure(varargin)
% CONFIGURE Installation of Mex files
if nargin == 0
    conf_Polylib = true;
    conf_Limit1d = true;
    conf_Limit2d = true;
    conf_SWE2d   = true;
    conf_Mesh    = true;
    conf_Post    = true;
else
    conf_Polylib = false;
    conf_Limit1d = false;
    conf_Limit2d = false;
    conf_SWE2d   = false;
    conf_Mesh    = false;
    conf_Post    = false;
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
            case 'Mesh'
                conf_Mesh = true;
            case 'Post'
                conf_Post = true;
            otherwise
                error(['No options for %s, please choose one of\n',...
                    '  Polylib\n  Limiter1d\n  Limiter2d\n  SWE2d\n',...
                    '  Mesh\n  Post\n'], varargin{i})
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

%% Mesh
if conf_Mesh
    localpath = '+Utilities/+Mesh/Mex';
    installpath = '+Utilities/+Mesh/';
    src = {'ResortVertex_Mex.c'};
    libsrc = {'VertexSort.c'};
    fprintf('==============Mesh=================\n')
    install(localpath, installpath, src, libsrc);
    fprintf('==============Mesh=================\n\n')
end% if

%% Limiter1d
if conf_Limit1d
    localpath = '+Utilities/+Limiter/Mex';
    installpath = '+Utilities/+Limiter/+Limiter1D';
    src = {'Minmod1d_Mex.c'};
    libsrc = {'Limiter.c'};
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
        'TVB_tri2d_Mex.c', ....
        'TVB_quad2d_Mex.c',...
        'VB2d_VA_Mex.c',...
        'VB2d_JK_Mex.c',...
        'VB2d_HWENO_Mex.c',...
        'TVB_detector2d_Mex.c'};
    libsrc = {'Limiter.c'};
    fprintf('==============Limiter2d=================\n')
    install(localpath, installpath, src, libsrc);
    fprintf('==============Limiter2d=================\n\n')
end
%% SWE2d
if conf_SWE2d
    localpath = 'SWE2d/Mex';
    installpath = localpath;
    src = {'SWE_Mex_BC2d.c', ...
        'SWE_Mex_Flux2d.c',...
        'SWE_Mex_HLLFlux2d.c',...
        'SWE_Mex_PositivePreserving2d.c'};
    libsrc = {'SWE2d.c'};
    fprintf('==============SWE2d=================\n')
    install(localpath, installpath, src, libsrc);
    fprintf('==============SWE2d=================\n\n')
end

%% Postprocess
if conf_Post
    localpath = '+Utilities/+PostProcess/Mex';
    installpath = '+Utilities/+PostProcess';
    src = {'FindLocCell_Mex.c'};
    libsrc = {'VectorOperator.c'};
    fprintf('==============Postprocess=================\n')
    install(localpath, installpath, src, libsrc);
    fprintf('==============Postprocess=================\n\n')
end% if
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
    mex('CFLAGS="$CFLAGS -Wall"','-O', file{:}, '-outdir', dirPath);
end% for

cd(pwdPath);
end