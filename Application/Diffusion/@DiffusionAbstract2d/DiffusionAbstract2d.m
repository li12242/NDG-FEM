classdef DiffusionAbstract2d < NdgPhysMat
    %DIFFUSIONABSTRACT2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties ( Constant )
        %> number of physical field [ C, tau1, tau2 ]
        Nfield = 3
        %> number of variable field
        Nvar = 1;
        %> index of variable in physical field
        varFieldIndex = 1;
    end
    
    properties ( Abstract, Constant )
        % constant viscosity
        miu
    end
    
    properties ( SetAccess = protected )
        % num of elements
        M
        % maximum order
        N
    end
    
    methods
        function obj = DiffusionAbstract2d()
            obj = obj@NdgPhysMat();
        end
    end
    
    methods( Hidden )
        
        function initPhysFromOptions( obj, mesh )
            initPhysFromOptions@NdgPhysMat( obj, mesh );
            finalTime = obj.getOption('finalTime');
            % for m = 1:obj.Nmesh
            %    mesh = obj.meshUnion(m);
            %    obj.fext{m} = obj.getExtFunc(mesh.x, mesh.y, finalTime);
            % end
            
            obj.viscositySolver = NdgQuadFreeStrongCentralVisSolver2d( obj, 1, 1 );
            for m = 1:obj.Nmesh % constant viscosity
                obj.viscositySolver.mx{1} = obj.miu;
                obj.viscositySolver.my{1} = obj.miu;
            end
        end
        
    end
    
    methods( Access = protected )
        function matEvaluateRHS( obj, fphys )
            obj.frhs{1}(:) = 0;
            obj.viscositySolver.matEvaluateRHS( fphys );
        end
    end
    
end

