function configure(varargin)
%CONFIGURE compile all the mex files.
%   Multiple mex files are used in the NDG-FEM model to speedup the
%   program, the user has to first compile these mex files and then use the
%   model to compute various applications.
%   Many mex files implement the OpenMP library to enable multi-threads
%   computating, the user can set the number of threads based on their own
%   computer devices as
% 
%       configure(4); % use 4 threads
%
%   The user also can use their own compiler 
%   
%   To clear all the compiled mex files, use follow command
%
%       configure -clear;
%
%   Author: li12242, Tianjin University.

delete_model = false;
switch nargin
    case 1
        compiler_path = '';
        
        if( isa(varargin{1}, 'double') )
            with_omp = varargin{1};
        elseif( isa(varargin{1}, 'char') )
            with_omp = 0;
            
            if( strcmp(varargin{1}, '-clear') )
                delete_model = true;
            end
        end
    case 2
        compiler_path = varargin{1};
        with_omp = varargin{2};
end

[compiler, cflags, ldflags] = set_the_compiler(compiler_path, with_omp);
%% +Polylib
path = '+Polylib';
srcfile = {'zwglj.c', 'zwgl.c', 'JacobiP.c', 'jacobfd.c', ...
    'GradJacobiP.c', 'Dglj.c'};
libfile = {'polylib.c'};
file_operater(path, srcfile, libfile, compiler, cflags, ldflags, delete_model);

%% +ndg_lib/+mesh
path = '+ndg_lib/+mesh/@mesh/private';
srcfile = {'cell_mean.c'};
libfile = {};
file_operater(path, srcfile, libfile, compiler, cflags, ldflags, delete_model);

path = '+ndg_lib/+mesh/@mesh2d/private';
srcfile = {'resort_vert.c'};
libfile = {'vert_sort.c'};
file_operater(path, srcfile, libfile, compiler, cflags, ldflags, delete_model);

%% +ndg_utility/+detector
path = '+ndg_utility/+detector/@detector2d/private';
srcfile = {'find_loc_cell.c'};
libfile = {'vec_operator.c'};
file_operater(path, srcfile, libfile, compiler, cflags, ldflags, delete_model);

%% +ndg_utility/+limiter
path = '+ndg_utility/+limiter/+VB/@VB_2d/private';
srcfile = {'vb_weno.c', 'vb_va.c', 'vb_jk.c', 'vertex_average.c', 'vertex_extreme.c'};
libfile = {'vb.c'};
file_operater(path, srcfile, libfile, compiler, cflags, ldflags, delete_model);

path = '+ndg_utility/+limiter/+BJ/private';
srcfile = {'vertex_extreme.c'};
libfile = {};
file_operater(path, srcfile, libfile, compiler, cflags, ldflags, delete_model);

path = '+ndg_utility/+limiter/+BJ/@BJ_2d/private';
srcfile = {'BJ_limit_2d.c'};
libfile = {};
file_operater(path, srcfile, libfile, compiler, cflags, ldflags, delete_model);

%% Conv2d
path = 'Conv2d/@conv2d/private';
srcfile = {'upwind_flux.c'};
libfile = {'conv2d.c'};
file_operater(path, srcfile, libfile, compiler, cflags, ldflags, delete_model);
path = 'Conv2d/@conv2d_adv_gq/private';
file_operater(path, srcfile, libfile, compiler, cflags, ldflags, delete_model);

path = 'Conv2d/@conv2d_refine_fv/private';
srcfile = {'rhs_term.c'};
libfile = {'conv2d.c', 'flux_term.c', 'surf_term.c'};
file_operater(path, srcfile, libfile, compiler, cflags, ldflags, delete_model);

srcfile = {'rhs_fv_term.c'};
libfile = {'conv2d.c'};
file_operater(path, srcfile, libfile, compiler, cflags, ldflags, delete_model);

%% NConv2d
path = 'NConv2d/@nconv2d/private';
srcfile = {'lf_flux.c'};
libfile = {'nconv2d.c'};
file_operater(path, srcfile, libfile, compiler, cflags, ldflags, delete_model);

%% SWE2d
path = 'SWE2d/@swe2d/private';
srcfile = {'hll_flux.c', 'nodal_flux.c', 'ppreserve.c', 'lf_flux.c'};
libfile = {'swe.c', 'bound_cond.c'};
file_operater(path, srcfile, libfile, compiler, cflags, ldflags, delete_model);

end% func

function file_operater(path, srcfile, libfile, ...
    compiler, cflags, ldflags, ...
    delete_model)

if delete_model
    clear_model(path, srcfile);
else
    install(path, srcfile, libfile, compiler, cflags, ldflags);
end
end% func

function [compiler, cflags, ldflags] = set_the_compiler(compiler_path, with_omp)
switch computer('arch')
    case 'maci64'
        if isempty(compiler_path) % default compiler path
            compiler = 'CC=/opt/intel/composer_xe_2015.5.222/bin/intel64/icc';
        else
            compiler = ['CC=',compiler_path, '/icc'];
        end
        
        if with_omp
            cflags = ['CFLAGS=$CFLAGS, -std=c99 -Wall -qopenmp -DDG_THREADS=', ...
                num2str(with_omp)];
            ldflags = {'-liomp5', '-largeArrayDims', '-lmwblas'};
        else
            cflags = 'CFLAGS=$CFLAGS, -std=c99 -Wall';
            ldflags = {'-largeArrayDims', '-lmwblas'};
        end
    case 'win64'
        compiler = '';
        if with_omp
            cflags = ['CFLAGS=$CFLAGS -fopenmp -DDG_THREADS=',...
                num2str(with_omp)];
            ldflags = {'LDFLAGS=$LDFLAGS -fopenmp', '-largeArrayDims', ...
                '-lmwblas'};
        else
            cflags = 'CFLAGS=$CFLAGS, -Wall';
            ldflags = {'-largeArrayDims', '-lmwblas'};
        end
    case 'glnxa64'
    otherwise
        error('configure error: Unknown platform.');
end
end% func

function clear_model(path, src)
fullPath = fullfile(pwd, path);
srcfile  = fullfile(fullPath, src);
fold = dir(fullPath);
obj_exten = mex_extern();
for n = 1:numel(srcfile)
    [ ~, name, ~ ] = fileparts( srcfile{n} );
    obj_file = [name, obj_exten];
    
    for i = 1:numel(fold)
        if ( strcmp(fold(i).name, obj_file) ) % find the compiled mex file
            obj_file = fullfile(fullPath, obj_file);
            delete(obj_file);
            break;
        end
    end
end% for

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
ndg_utility.cprintf('string', ...
    ['%s\nCFLAGS=%s\nLDFLAGS=', sformat, '\n'], ...
    compiler, cflags, ldflags{:});

for i = 1:numel(srcfile)
    if ( iscompiled(srcfile{i}, libfile) ) 
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

function tf = iscompiled( src_file, lib_files )
% Check whether the source files and library files are changed.
tf = false;
obj_exten = mex_extern();
[ path,name,src_exten ] = fileparts( src_file );
src_file_name = [ name, src_exten ];
obj_file_name = [ name, obj_exten ];
fold = dir(path);
obj_datenum = 0;
src_datenum = 0;
lib_datenum = 0;
for i = 1:numel( fold )
    if strcmp( fold(i).name, obj_file_name )
        obj_datenum = fold(i).datenum;
    elseif strcmp( fold(i).name, src_file_name )
        src_datenum = fold(i).datenum;
    end
    
    for j = 1:numel(lib_files)
        [ ~,name,src_exten ] = fileparts( lib_files{j} );
        lib_file_name = [ name,src_exten ];
        if strcmp( fold(i).name, lib_file_name )
            lib_datenum = max(lib_datenum, fold(i).datenum);
        end
    end
end% for

if obj_datenum % the obj file exists
    if (src_datenum < obj_datenum) && (lib_datenum < obj_datenum)
        % the object file is new
        tf = true;
    end% 
end% if

end% func

function obj_exten = mex_extern()
switch computer('arch')
    case 'maci64'
        obj_exten = '.mexmaci64';
    case 'win64'
        obj_exten = '.mexw64';
    case 'glnxa64'
        obj_exten = '.mexa64';
end
end% func

