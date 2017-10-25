classdef Adv2d < NdgPhys
    %ADV2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        %> number of variable field
        Nvar = 1
        %> number of physical field
        Nfield = 3
        %> the variable field
        varFieldId = 1
    end
    
    methods
        function obj = Adv2d(solver, mesh)
            obj = obj@NdgPhys(solver, mesh);
        end
        
        function mxSolve(obj)
            obj.initPhysField;
            obj.fphy = AdvSolver2d(obj);
        end
        
        function drawResult(obj, varargin)
            if( nargin == 1 )
                deltaTimeStep = 1;
            else
                deltaTimeStep = varargin{1};
            end
            file = [obj.solver.outputNetcdfCaseName, '.0-1.nc'];
            time = ncread(file, 'time');
            Ntime = numel(time);

            for t = 1:deltaTimeStep:Ntime
                field = ncread(file, 'fphy', [1,1,1,t], [inf, inf, 1, 1]);
                obj.mesh.draw(field);
                zlim([-0.2, 1.2]);
                drawnow;
            end
        end
    end
    
end

