classdef SWEAbstract2d < NdgPhysMat
    
    properties(Abstract, Constant)
        %> wet/dry depth threshold
        hmin
        %> gravity acceleration
        gra
        %> number of physical field
        Nfield
        %> number of variable field
        Nvar
        %> index of variable in physical field
        varFieldIndex
    end
    
    properties
        %> gradient of bottom elevation
        zGrad
        %>
        coriolisSolver
        %>
        frictionSolver
        %>
        windSolver
        %>
        wind
    end
    
    methods
        draw( obj, varargin )
    end
    
    methods( Hidden )
        
        function [ fM, fP ] = matEvaluateSurfaceValue( obj, mesh, fphys, fext )
            [ fM, fP ] = mxEvaluateSurfaceValue( obj.hmin, obj.gra, ...
                mesh.eidM, mesh.eidP, mesh.nx, mesh.ny, mesh.eidtype, ...
                fphys, fext, mesh.EToE, mesh.EToR );
        end
        
        function [ fluxM ] = matEvaluateSurfFlux( obj, mesh, nx, ny, fm )
            [ fluxM ] = mxEvaluateSurfFlux( obj.hmin, obj.gra, nx, ny, fm);
        end
        
        function [ fluxS ] = matEvaluateSurfNumFlux( obj, mesh, nx, ny, fm, fp )
            [ fluxS ] = matEvaluateHLLNumFlux2d( obj, nx, ny, fm, fp );
        end
    end
    
    methods( Abstract, Access = protected )
        %> 
        matEvaluateTopographySourceTerm( obj, fphys )
    end
    
    methods( Access = protected )
        %>
        function [ fphys ] = matEvaluatePostFunc(obj, fphys)
            for m = 1:obj.Nmesh
                hc = obj.meshUnion(m).GetMeshAverageValue( fphys{m}(:,:,1) );
                qxc = obj.meshUnion(m).GetMeshAverageValue( fphys{m}(:,:,2) );
                qyc = obj.meshUnion(m).GetMeshAverageValue( fphys{m}(:,:,3) );
                fphys{m}(:,:,1:3) = mxEvaluatePostFunc2d( obj.hmin, fphys{m}, hc, qxc, qyc );
            end
            obj.matUpdateWetDryState( fphys );
        end
        
        function matUpdateWetDryState(obj, fphys)
            for m = 1:obj.Nmesh
                wetflag = all( fphys{m}(:,:,1) > obj.hmin );
                obj.meshUnion(m).EToR( ~wetflag ) = int8( NdgRegionType.Dry );
                obj.meshUnion(m).EToR(  wetflag ) = int8( NdgRegionType.Wet );
            end
        end
        
        function [ dt ] = matUpdateTimeInterval( obj, fphys )
            dt = nan;
            for m = 1:obj.Nmesh
                dx = sqrt( obj.meshUnion(m).LAV );
                N = obj.meshUnion(m).cell.N;
                dtm = mxUpdateTimeInterval2d( ...
                    obj.hmin, ...
                    obj.gra, ...
                    N, ...
                    dx, ...
                    obj.meshUnion(m).EToR, ...
                    fphys{m} );
                
                if ( dtm > 0 )
                    dt = min(dt, dtm * obj.cfl);
                end
            end
        end
        
        function matEvaluateSourceTerm( obj, fphys )
            % frhs = frhs + BottomTerm
            obj.matEvaluateTopographySourceTerm( fphys );
            
            % frhs = frhs + CoriolisTerm
            obj.coriolisSolver.evaluateCoriolisTermRHS(obj, fphys);
            
            % frhs = frhs + FrictionTerm
            obj.frictionSolver.evaluateFrictionTermRHS(obj, fphys);
            
            % frhs = frhs + WindTerm
            obj.windSolver.evaluateWindTermRHS(obj, fphys);
        end
        
        function [ fluxS ] = matEvaluateHLLNumFlux2d( obj, nx, ny, fm, fp )
            [ fluxS ] = mxEvaluateHLLNumFlux2d( obj.hmin, obj.gra, nx, ny, fm, fp );
        end
        
        function [ fluxS ] = matEvaluateLFNumFlux2d( obj, mesh, nx, ny, fm, fp )
            
        end
    end
    
end
