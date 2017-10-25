classdef NdgPhys < handle
    %NDGPHYS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(SetAccess = protected )
        %> Number of mesh 
        Nmesh
        %> mesh objects
        mesh
        %> solver setting
        solver
    end
    
    properties(Constant, Abstract)
        %> number of variable field
        Nvar
        %> number of physical field
        Nfield
        %> final time
        finalTime
        %> start time
        initialTime
        %> index of variable field
        varFieldId
    end
    
    properties
        %> field values, dimensions - [Np, K, Nfield, Nmesh]
        fphy
    end
    
    methods
        function obj = NdgPhys(solver, mesh)
            obj.solver = solver;
            obj.Nmesh = numel(mesh);
            obj.mesh = mesh;
        end
    end
    
end

