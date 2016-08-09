classdef SpongeBC1d < MultiRegions.BoundCondition.SpongeBC
    %SPONGEBC1D Summary of this class goes here
    %   Detailed explanation goes here
    
properties
    SpNToBVCoeff % coefficient for all nodes in sponge elements to BV
end

methods
    function obj = SpongeBC1d(mesh, BCflag, fileName)
        obj = obj@MultiRegions.BoundCondition.SpongeBC(BCflag, fileName);
        % the connection between sponge element to the boundary vertex
        obj.SpEToBV = zeros(obj.nSpE, 1);
        xbv = obj.BCfile.GetVarData('xb');
        xp  = mesh.x(:, obj.SpEToE);
        xpc = mean(xp); % centre coordinate of sponge element
    end% func
end
    
end

