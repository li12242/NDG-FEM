classdef SpongeBC
%SONGEBC Sponge boundary condition.
%   Detailed explanation goes here

properties
    BCfile          % boundary condition files
    time            % time vector
    nSpE            % number of sponge element
    SpEToE          % sponge element to computation elements
end

methods
    %% SpongeBC 
    % construction function of the sponge boundary condition
    function obj = SpongeBC(BCflag, fileName)
        % boundary netcdf file
        obj.BCfile = Utilities.PostProcess.ResultFile(fileName);
        % get nBV
        ncid  = netcdf.open(fileName, 'NC_NOWRITE');
        dimid = netcdf.inqDimID(ncid, 'ne'); % the dimension name must be ne
        [~, ne] = netcdf.inqDim(ncid, dimid);
        % get time
        obj.time = obj.BCfile.GetVarData('time');
        % the connection between sponge element to computation element
        obj.nSpE = sum(BCflag);
        obj.SpEToE = find(BCflag);
        % check for number of sponge layer element
        if ne ~= obj.nSpE
            error(['The number of sponge layer elements in OBC file [',...
                num2str(ne),...
                '] does not equal to the input variable BCflag [', ...
                num2str(obj.nSpE) ,']']);
        end% if
    end% func
end
    
end

