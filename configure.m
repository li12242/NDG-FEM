%======================================================================
%> @brief Compile all the mex files.
%>
%> Multiple mex files are used in the NDG-FEM model to speedup the 
%> computation, the user has to first compile these mex files and then use 
%> the model to simulate various applications.
%> Many mex files implement the OpenMP library to enable multi-threads
%> computating, the user can set the number of threads based on their own
%> computer devices as
%>
%> @code 
%>   configure(4); % use 4 threads
%> @endcode
%> 
%> The user also can use their own compiler 
%>   
%> To remove all the compiled mex files, use following command:
%> @code
%>   configure -clear;
%> @endcode
%> 
%======================================================================
%> This class is part of the NDG-FEM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
%======================================================================
function configure(varargin)

deleteMod = false;
compilerPath = '';

if (nargin == 1) && ( isa(varargin{1}, 'double') )
    threadNum = varargin{1};
elseif (nargin == 1) && ( isa(varargin{1}, 'char') )
    threadNum = 0;
    compilerPath = varargin{1};
    if ( strcmp(varargin{1}, '-clear') )
        deleteMod = true;
    end
elseif (nargin == 2)
    compilerPath = varargin{1};
    threadNum = varargin{2};
end

[compiler, cflags, ldflags] = set_the_compiler( compilerPath, threadNum );

% Polylib
path = 'thirdParty/Polylib/';
srcfile = {[path, 'zwglj.c'], ...
    [path, 'zwgl.c'], ... 
    [path, 'JacobiP.c'], ...
    [path, 'GradJacobiP.c'] };
libfile = { [path, 'polylib.c'] };
outPath = 'lib/';
file_operater(outPath, srcfile, libfile, compiler, cflags, ldflags, deleteMod);

% Advection
srcfile = {'Advection/@AdvUniformUnion2d/src/AdvSolver2d.c'};
libfile = {'NdgCell/src/NdgCell.c', ...
    'NdgCell/src/NdgCellIO.c', ...
    'NdgMesh/src/NdgMesh.c', ... 
    'NdgMesh/src/NdgMeshIO.c', ...
    'NdgPhys/src/NdgPhys.c', ...
    'NdgPhys/src/NdgPhysIO.c', ...
    'NdgPhys/src/NdgPhysSolver.c', ...
    'NdgPhys/src/NdgTerporalDiscrete.c', ...
    'NdgPhys/src/mxNdgPhys.c', ...
    'NdgUtils/src/MatUtils.c',...
    'NdgUtils/src/NdgNetCDF.c',...
    'NdgUtils/src/Utils.c',...
    'NdgUtils/src/NdgLimiter.c'};
outPath = 'lib/';
cflags = [cflags, ' -DPROFILE ', ...
    ' -INdgCell/src/ -INdgMesh/src/ -INdgPhys/src -INdgUtils/src '];
ldflags = [ldflags, ' -lnetcdf '];
file_operater(outPath, srcfile, libfile, compiler, cflags, ldflags, deleteMod);

end% func

function file_operater(outpath, srcfile, libfile, ...
    compiler, cflags, ldflags, deleteModel)
if deleteModel
    remove_lib_file(outpath, srcfile);
else
    install_lib_file(outpath, srcfile, libfile, compiler, cflags, ldflags);
end
end% func

%> set the compile options.
function [compiler, cflags, ldflags] ...
    = set_the_compiler( compilerPath, Nthread )

compiler = '';
cflags = 'CFLAGS=$CFLAGS -std=c99 -Wall ';
ldflags = 'LDFLAGS=$LDFLAGS -lmwblas ';

switch computer('arch')
    case 'maci64'
        if isempty(compilerPath) % default compiler path
            compiler = ...
                'CC=/opt/intel/composer_xe_2015.5.222/bin/intel64/icc';
        else
            compiler = ['CC=',compilerPath, '/icc'];
        end
        
        if Nthread
            cflags = [cflags, ' -qopenmp -DDG_THREADS=', num2str(Nthread)];
            ldflags = [ldflags, '-liomp5 '];
        end
    case 'win64'
        if Nthread
            cflags = [cflags,' -fopenmp -DDG_THREADS=',num2str(Nthread)];
            ldflags = [ldflags, ' -fopenmp'];
        end
    case 'glnxa64'
    otherwise
        msgID = 'configure:set_the_compiler';
        msgtext = 'Unknown computer architecture.';
        ME = MException(msgID, msgtext);
        throw(ME);
end
end% func

%> remove the compiled object files.
function remove_lib_file( path, srcfile )
libPath = [pwd, '/', path];
ext = mex_extern();
fold = dir(libPath);
for n = 1:numel(srcfile)
    [ ~, name, ~ ] = fileparts( srcfile{n} );
    obj_file = [name, ext];
    for i = 1:numel(fold)
        if ( strcmp(fold(i).name, obj_file) ) % find the compiled mex file
            delete( fullfile(libPath, obj_file) );
            break;
        end
    end
end% for

end% func

%> compile the source file and install the object files.
function install_lib_file(outpath, src, libsrc, compiler, cflags, ldflags)
outDir = [pwd, '/', outpath];
fprintf('installing source files to %s.\n', outpath);
fprintf('%s\nCFLAGS=%s\nLDFLAGS=%s\n', compiler, cflags, ldflags);

for i = 1:numel(src)
    if ( isCompiled(outDir, src{i}, libsrc) ) 
        continue; 
    end
    fprintf('\n%s...\n', src{i});
    file = [src(i), libsrc{:}];
    mex(compiler, cflags, '-O', '-largeArrayDims', ...
        ldflags, file{:}, '-outdir', outDir);
end% for
fprintf('finish installing srouce files.\n');
end

%> check whether the source and library files are changed.
function tf = isCompiled( outPath, srcFile, libSrdFile )
tf = false;
obj_exten = mex_extern();
[ ~,name,src_exten ] = fileparts( srcFile );
src_file_name = [ name, src_exten ];
obj_file_name = [ name, obj_exten ];
fold = dir( outPath );
obj_datenum = 0;
src_datenum = 0;
lib_datenum = 0;
for i = 1:numel( fold )
    if strcmp( fold(i).name, obj_file_name )
        obj_datenum = fold(i).datenum;
    elseif strcmp( fold(i).name, src_file_name )
        src_datenum = fold(i).datenum;
    end
    
    for j = 1:numel(libSrdFile)
        [ ~,name,src_exten ] = fileparts( libSrdFile{j} );
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

%> return the file suffixion based on different platform
function ext = mex_extern()
switch computer('arch')
    case 'maci64'
        ext = '.mexmaci64';
    case 'win64'
        ext = '.mexw64';
    case 'glnxa64'
        ext = '.mexa64';
end
end% func

