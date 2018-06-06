classdef NdgNcDim < handle
    
    properties(SetAccess=private)
        %> dimension name
        name
        %> positive integer, 0 for unlimited
        length
        %> dimension id in NetCDF file
        id
    end
    
    methods
        function [ obj ] = NdgNcDim( name, length )
            obj.name = name;
            obj.length = length;
        end% func
        
        function setDimId( obj, dimId )
            [ obj.id ] = dimId;
        end
        
        function defineIntoNetcdfFile(obj, ncid)
            if (obj.length > 0)
                obj.id = netcdf.defDim(ncid, obj.name, obj.length);
            else
                typeid = netcdf.getConstant( 'NC_UNLIMITED' );
                obj.id = netcdf.defDim(ncid, obj.name, typeid);
            end
        end% func
        
    end% methods
    
end

