classdef NdgNcFile < handle
    
    properties(SetAccess = protected)
        %> true for file is open
        isOpen = false
        %> ID of ncfile
        ncid
        %> array of dimensions in NetCDF file
        ncDim
        %> array of variables in NetCDF file
        ncVar
        %> file name of NetCDF file
        fileName
    end
    
    methods
        %======================================================================
        %> @brief Brief description of the function
        %>
        %> More detailed description.
        %>
        %> @param arg1 First argument
        %> @param arg2 Second argument
        %>
        %> @retval out1 return value for the first output variable
        %> @retval out2 return value for the second output variable
        %======================================================================
        %> This function is part of the NDGOM software.
        %> @author li12242, Tianjin University, li12242@tju.edu.cn
        %======================================================================
        function obj = NdgNcFile( filename, ncdim, ncvar )
            obj.fileName = filename;
            obj.ncDim = ncdim;
            obj.ncVar = ncvar;
            obj.isOpen = false;
        end% func
        
        function delete( obj )
            if(obj.isOpen) % if netcdf file is still open
                obj.isOpen = false;
                netcdf.close( obj.ncid );
            end
        end% func
        
        function defineIntoNetcdfFile( obj )
            obj.ncid = netcdf.create( obj.fileName, 'CLOBBER');
            obj.isOpen = true;
            for n = 1:numel(obj.ncDim)
                obj.ncDim(n).defineIntoNetcdfFile( obj.ncid );
            end
            
            for n = 1:numel(obj.ncVar)
                obj.ncVar(n).defineIntoNetcdfFile( obj.ncid );
            end
            netcdf.endDef(obj.ncid);
        end
    end
    
end

function [ncid, ncDim, ncVar] = readFromNetcdfFile( fileName )
ncid = netcdf.open(fileName, 'NOWRITE');

dimId = netcdf.inqDimIDs( ncid );
Ndim = numel(dimId);
ncdim = cell(Ndim, 1);
for n = 1:Ndim
    % read dimension name and ID
    [ dimName, dimLength ] = netcdf.inqDim( ncid, dimId(n) );
    dim = NdgNcDim( dimName, dimLength );
    dim.setDimId( dimId(n) );
    ncdim{n} = dim;
end
ncDim = [ncdim{:}];

varId = netcdf.inqVarIDs( ncid );
Nvar = numel(varId);
ncvar = cell(Nvar, 1);
for n = 1:Nvar
    % read variable name and ID
    [varName, dataType, dimIds, ~] = netcdf.inqVar(ncid, varId(n) );
    ncdim = ncDim( dimIds + 1 );
    var = NdgNcVar( varName, ncdim, NdgNcType( dataType ) );
    var.setVarId( varId(n) );
    ncvar{n} = var;
end
ncVar = [ ncvar{:} ];

end% func
