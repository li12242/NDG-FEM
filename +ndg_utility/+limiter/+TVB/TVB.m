classdef TVB
    %TVB Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(SetAccess = protected)
        mesh
        cell
        V       % Vandermonde matrix
        invV    % inverse of Vandermonde matrix
        h2      % square of each elements' length
        xc      % 
        yc
        zc
    end
    
    methods
        function obj = TVB(mesh, cell)
            obj.mesh = mesh;
            obj.cell = cell;
            obj.V = cell.V;
            obj.invV = inv(obj.V);
            
            obj.h2 = obj.mesh.vol.^2; % h2
            obj.xc = obj.mesh.cell_mean(obj.mesh.x);
            obj.yc = obj.mesh.cell_mean(obj.mesh.y);
            obj.zc = obj.mesh.cell_mean(obj.mesh.z);
        end
        
        
    end
    
end

