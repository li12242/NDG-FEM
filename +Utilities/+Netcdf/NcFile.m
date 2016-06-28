classdef NcFile < handle
    properties
        dimArray % array of dimensions
        varArray % array of variables
        name % file name
        fileID % file id
    end% properties
    
    methods
        %% NcFile
        % Define NetCDF output file structure.
        %
        % Usages:
        %   time = Utilities.netcdf.NcDim('time', 0); % unlimited dimensions
        %   node = Utilities.netcdf.NcDim('node', mesh.nNode);
        %
        %   x    = Utilities.netcdf.NcVar('x',    node, 'double');
        %   y    = Utilities.netcdf.NcVar('y',    node, 'double');
        %   t    = Utilities.netcdf.NcVar('time', time, 'double');
        %   var  = Utilities.netcdf.NcVar('var',  [node, time], 'double');
        %   file = Utilities.netcdf.NcFile('file', [node, time], [x, y, t, var]);
        %
        function obj = NcFile(name, dimArray, varArray)
            obj.name = name; % no suffix 
            obj.dimArray = dimArray;
            obj.varArray = varArray;
        end% func
        
        %% CreateFile
        % Generate and open NetCDF files based on object.
        %
        % Attention:
        %
        %   Remember to close the file with CloseFile
        %
        % Usages:
        %   file.CreateFile;
        function CreateFile(obj)
            fileName = [obj.name, '.nc'];
            obj.fileID = netcdf.create(fileName,'CLOBBER');

            % Define dimensions
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
            
            % Define variables
            nvar = numel(obj.varArray);
            
            for i = 1:nvar
                obj.varArray(i).id = netcdf.defVar(obj.fileID, ...
                    obj.varArray(i).name, ...
                    obj.varArray(i).type, ...
                    [obj.varArray(i).dimArray(1:end).id]);
            end% for
            % end defination
            netcdf.endDef(obj.fileID);
        end% func
        
        %% PutVarPart
        % Put values into NetCDF files
        % Attention: the index of netcdf start form 0 (C style)
        function putVarPart(obj, varName, startIndex, len, varValue)        
            % get variable ID
            varID = getVarID(obj, varName);
            % put values
            netcdf.putVar(obj.fileID, varID, startIndex, len, varValue);
        end% func
        %% getVarID
        % Get variable ID based on its name
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
        
        %% CloseFile
        % Close NetCDF file
        function CloseFile(obj)
            netcdf.close(obj.fileID);
        end% func
        
    end% methods
end% class