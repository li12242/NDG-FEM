classdef NdgNcVar < handle

    properties(SetAccess=private)
        %> NdgNcDim object
        dims
        %> variable name
        name        
        %> variable ID
        id          
        %> variable types
        dataType
    end% properties
    
    methods
        function obj = NdgNcVar( name, ncDims, dataType )
            obj.name = name;
            obj.dims = checkNetcdfDims( ncDims );
            obj.dataType = checkDataType(dataType);
        end% func
        
        function setVarId( obj, varId )
            obj.id = varId;
        end% func
        
        function defineIntoNetcdfFile(obj, ncid)
            obj.id = netcdf.defVar(ncid, obj.name, obj.dataType, ...
                [obj.dims(1:end).id] );
        end% func
    end
    
end

function ncDim = checkNetcdfDims(ncDim)
% set dimension arrays
Ndim = numel(ncDim);

for n = 1:Ndim
    if ( n ~= Ndim ) && ( ncDim(n).length <= 0 )
        msgID = [mfilename, ':DimensionError.'];
        msgtext = ['The length of NdgNetcdfDim ', ...
            ncDim(n).name, ' should be positive.'];
        throw( MException(msgID, msgtext) );
    elseif ( ncDim(n).length < 0 )
        % the length of last dim can equal to zero
        msgID = [mfilename, ':DimensionError.'];
        msgtext = ['The length of last NdgNetcdfDim ', ...
            ncDim(n).name, ' should be positive or zero (Unlimited).'];
        throw( MException(msgID, msgtext) );
    end
end

end

function dataType = checkDataType( type )
% set variable types as double, float, short or int.
switch type
    case enumNcData.NC_DOUBLE
        dataType = 'double';
    case enumNcData.NC_FLOAT
        dataType = 'single';
    case enumNcData.NC_INT
        dataType = 'int';
    case enumNcData.NC_SHORT
        dataType = 'int16';
    case enumNcData.NC_CHAR
        dataType = 'char';
    case enumNcData.NC_BYTE
        dataType = 'int8';
    otherwise
        msgID = [mfilename, ':UnknownNcDataType.'];
        msgtext = 'The data type should be an enumeration of NdgNcDataType class.';
        throw( MException(msgID, msgtext) );
end% switch
end% func

