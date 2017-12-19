%> @brief Shallow water equations with solid continuous topography
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDGOM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef SWEAbstractCB2d < NdgPhysMat

    properties(Abstract, Constant)
        %> wet/dry depth threshold
        hmin
        %> gravity acceleration
        gra
        %> Physical field
        Nfield
    end
    
    properties(Constant)
        %> Variable field - {h, hu, hv}
        Nvar = 3
        %> field index of variable field
        varFieldIndex = [1,2,3]
    end
    
    properties
        %> RHS term Solvers
        coriolisSolver
        frictionSolver
    end
    
    methods
        function obj = SWEAbstractCB2d()
            obj = obj@NdgPhysMat();
        end
        
        drawSurfaceBot( obj )
        initPhysFromOptions( obj, mesh );
    end
    
    methods( Hidden )
        function [ E, G ] = matEvaluateFlux( obj, mesh, fphys )
            [ E, G ] = mxEvaluateFlux2d( obj.hmin, obj.gra, mesh.EToR, fphys );
        end
        
        function [ fluxS ] = matEvaluateSurfNumFlux( obj, mesh, nx, ny, fphys, fext )
            [ fluxS ] = mxEvaluateSurfNumFlux( obj.hmin, obj.gra, ...
                mesh.eidM, mesh.eidP, mesh.eidtype, nx, ny, fext, fphys );
        end
        
        function [ flux ] = matEvaluateSurfFlux( obj, mesh, nx, ny, fphys )
            [ flux ] = mxEvaluateSurfFlux( obj.hmin, obj.gra, ...
                mesh.eidM, mesh.eidtype, nx, ny, fphys);
        end
    end
    
    methods( Access = protected )
         
        [ fphys ] = matEvaluatePostFunc(obj, fphys)
        [ dt ] = matUpdateTimeInterval( obj, fphys )
        matEvaluateSourceTerm( obj, mesh, fphys )
        
        function matUpdateWetDryState(obj, fphys)
            for m = 1:obj.Nmesh
                wetflag = all( fphys{m}(:,:,1) > obj.hmin ); 
                obj.meshUnion(m).EToR( ~wetflag ) = int8( NdgRegionType.Dry );
                obj.meshUnion(m).EToR(  wetflag ) = int8( NdgRegionType.Wet );
            end
        end
    end
    
end

