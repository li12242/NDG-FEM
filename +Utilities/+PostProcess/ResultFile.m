classdef ResultFile < handle
    properties
        filename
    end% properties
    
    methods
        function obj = ResultFile(filename)
            obj.filename = filename;
        end% func
        
        %% Get variable data 
        function data = GetVarData(obj, varName)
            data = ncread(obj.filename, varName);
        end% func
        
        %% Get variable data at spicific time step
        function data = GetTimeVarData(obj, varName, ist)
            % get number of dimensions
            ncid   = netcdf.open(obj.filename);
            varid  = netcdf.inqVarID(ncid, varName);
            [~, ~, dimids, ~] = netcdf.inqVar(ncid, varid);
            ndim   = numel(dimids);
            % get each dimension length
            dimlen = zeros(ndim, 1);
            for i  = 1:ndim
                [~, dimlen(i)] = netcdf.inqDim(ncid, dimids(i));
            end
            % get result
            % warnning: the index is C style, start from 0, and the last
            % dimension must be time.
            start  = [zeros(ndim-1, 1); ist-1];
            count  = [dimlen(1:ndim-1); 1];
            data   = netcdf.getVar(ncid, varid, start, count);
            netcdf.close(ncid);
        end% func
        
    end% methods
end% classdef