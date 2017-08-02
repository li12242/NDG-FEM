function configure(varargin)
% CONFIGURE 编译 mex 文件
% 输入参数包括
%   1. 编译器路径（可省略）；
%   2. 是否选择 openmp 编译，0-不使用， 其他-OpenMP 线程数；
switch nargin
    case 1
        compiler_path = '';
        with_omp = varargin{1};
    case 2
        compiler_path = varargin{1};
        with_omp = varargin{2};
end
switch computer('arch')
    case 'maci64'
        if isempty(compiler_path) % default compiler path
            compiler = 'CC=/opt/intel/composer_xe_2015.5.222/bin/intel64/icc';
        else
            compiler = ['CC=',compiler_path, '/icc'];
        end
        
        if with_omp
            cflags = ['CFLAGS=$CFLAGS, -Wall -qopenmp -DDG_THREADS=', ...
                num2str(with_omp)];
            ldflags = {'-liomp5', '-largeArrayDims', '-lmwblas'};
        else
            cflags = 'CFLAGS=$CFLAGS, -Wall';
            ldflags = {'-largeArrayDims', '-lmwblas'};
        end
    case 'win64'
        compiler = '';
        if with_omp
            cflags = ['CFLAGS=$CFLAGS -fopenmp -DDG_THREADS=',...
                num2str(with_omp)];
            ldflags = {'-liomp5', '-largeArrayDims', '-lmwblas'};
        else
            cflags = 'CFLAGS=$CFLAGS, -Wall';
            ldflags = {'-largeArrayDims', '-lmwblas'};
        end
    case 'glnxa64'
    otherwise
        error('configure error: Unknown platform.');
end

%% +Polylib
path = '+Polylib';
srcfile = {'zwglj.c', 'zwgl.c', 'JacobiP.c', 'jacobfd.c', ...
    'GradJacobiP.c', 'Dglj.c'};
libfile = {'polylib.c'};
install(path, srcfile, libfile, compiler, cflags, ldflags);

%% +ndg_lib/+mesh
path = '+ndg_lib/+mesh/@mesh/private';
srcfile = {'cell_mean.c'};
libfile = {};
install(path, srcfile, libfile, compiler, cflags, ldflags);

path = '+ndg_lib/+mesh/@quad_mesh/private';
srcfile = {'resort_vert.c'};
libfile = {'VertexSort.c'};
install(path, srcfile, libfile, compiler, cflags, ldflags);
path = '+ndg_lib/+mesh/@tri_mesh/private';
install(path, srcfile, libfile, compiler, cflags, ldflags);

%% +ndg_utility/+detector
path = '+ndg_utility/+detector/@detector2d/private';
srcfile = {'find_loc_cell.c'};
libfile = {'vec_operator.c'};
install(path, srcfile, libfile, compiler, cflags, ldflags);

%% +ndg_utility/+limiter
path = '+ndg_utility/+limiter/+VB/@VB_2d/private';
srcfile = {'vb_weno.c', 'vertex_average.c', 'vertex_extreme.c'};
libfile = {'vb.c'};
install(path, srcfile, libfile, compiler, cflags, ldflags);

%% Conv2d
path = 'Conv2d/@conv2d/private';
srcfile = {'upwind_flux.c'};
libfile = {'conv2d.c', 'bc.c'};
install(path, srcfile, libfile, compiler, cflags, ldflags);
path = 'Conv2d/@conv2d_adv_gq/private';
install(path, srcfile, libfile, compiler, cflags, ldflags);

path = 'Conv2d/@conv2d_refine_fv/private';
srcfile = {'rhs_term.c'};
libfile = {'conv2d.c', 'flux_term.c', 'surf_term.c'};
install(path, srcfile, libfile, compiler, cflags, ldflags);

%% SWE2d
path = 'SWE2d/@swe2d/private';
srcfile = {'hll_flux.c', 'nodal_flux.c', 'ppreserve.c', 'lf_flux.c'};
libfile = {'swe.c', 'bound_cond.c'};
install(path, srcfile, libfile, compiler, cflags, ldflags);

end% func

%% install function

function install(path, src, libsrc, compiler, cflags, ldflags)
% Compile the source file and put into spicific directory.
addpath(pwd)
pwdPath  = pwd;
fullPath = fullfile(pwd, path);

srcfile  = fullfile(fullPath, src);
libfile  = fullfile(fullPath, libsrc);

cd(fullPath);

sformat = repmat('%s ', 1, numel(ldflags));
ndg_utility.cprintf('key', '=========installing %s=========\n', path);
ndg_utility.cprintf('string', ['%s\nCFLAGS=%s\nLDFLAGS=', sformat, '\n'], ...
    compiler, cflags, ldflags{:});

for i = 1:numel(srcfile)
    if ( iscompiled(srcfile) ) 
        continue; 
    end
    fprintf('\n%s/%s...\n', path,src{i});
    file = [srcfile(i), libfile{:}];
    mex(compiler, cflags, '-O', ldflags{:}, file{:});
end% for
ndg_utility.cprintf('key', ...
    '=========finish installing %s=========\n\n', path);
cd(pwdPath);
end

function tf = iscompiled(srfile)
tf = false;
end

