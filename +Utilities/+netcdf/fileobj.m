classdef fileobj < handle
    properties
        dimArray % array of dimensions
        varArray % array of variables
        name % file name
        fileID % file id
    end% properties
    
    methods
        function obj = fileobj(name, dimArray, varArray)
            obj.name = name;
            obj.dimArray = dimArray;
            obj.varArray = varArray;
        end% func
        
        function createFile(obj)
        % create file in Netcdf format
            fileName = [obj.name, '.nc'];
            obj.fileID = netcdf.create(fileName,'CLOBBER');

            % define dimensions
            ndim = numel(obj.dimArray);
            for i = 1:ndim
                if obj.dimArray(i).len > 0
                    obj.dimArray(i).id = netcdf.defDim(obj.fileID, ...
                        obj.dimArray(i).name, ...
                        obj.dimArray(i).len);
                elseif obj.dimArray(i).len <= 0
                    obj.dimArray(i).id = netcdf.defDim(obj.fileID, ...
                        obj.dimArray(i).name, ...
                        netcdf.getConstant('NC_UNLIMITED'));
                end% if
            end% for
            
            % define variables
            nvar = numel(obj.varArray);
            
            for i = 1:nvar
                obj.varArray(i).id = netcdf.defVar(obj.fileID, ...
                    obj.varArray(i).name, ...
                    obj.varArray(i).type, ...
                    [obj.varArray(i).dimArray(1:end).id]);
            end% for
            
            netcdf.endDef(obj.fileID);
        end% func
        
        function putVarPart(obj, varName, startIndex, len, varValue)
        % put values into netcdf files
        % Attention: the index of netcdf start form 0 (C style)
        
            % get variable ID
            varID = getVarID(obj, varName);
            
            % put values
            netcdf.putVar(obj.fileID, varID, startIndex, len, varValue);
        end% func
        
        function ID = getVarID(obj, varName)
        % get variable ID through its name
            nvar = numel(obj.varArray); % number of arrays
            for i = 1:nvar
                % compare variable names
                if strcmp(varName, obj.varArray(i).name)
                    ID = obj.varArray(i).id;
                    return;
                end% if
            end% for
        end% func
        
        function closeFile(obj)
            netcdf.close(obj.fileID);
        end% func
        
    end% methods
end% class