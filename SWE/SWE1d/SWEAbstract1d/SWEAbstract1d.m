classdef SWEAbstract1d < NdgPhysMat
    
    properties( Abstract, Constant )
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
        frictionSolver
    end
    
    methods( Hidden )
        
        function [ fM, fP ] = matEvaluateSurfaceValue( obj, mesh, fphys, fext )
            [ fM, fP ] = mxEvaluateSurfaceValue( obj.hmin, obj.gra, ...
                mesh.eidM, mesh.eidP, ...
                mesh.nx, mesh.eidtype, fphys, fext );
        end
        
        function [ fluxM ] = matEvaluateSurfFlux( obj, mesh, nx, fm )
            [ fluxM ] = mxEvaluateSurfFlux( obj.hmin, obj.gra, nx, fm);
        end
        
        function [ fluxS ] = matEvaluateSurfNumFlux( obj, mesh, nx, fm, fp )
            [ fluxS ] = matEvaluateHLLNumFlux1d( obj, nx, fm, fp );
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
                fphys{m}(:,:,1:2) = mxEvaluatePostFunc1d( obj.hmin, fphys{m}, hc, qxc );
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
                N = obj.meshUnion(m).cell.N;
                dtm = mxUpdateTimeInterval1d( obj.gra, N, obj.meshUnion(m).LAV, ...
                    obj.meshUnion(m).EToR, fphys{m} );
                if ( dtm > 0 )
                    dt = min(dt, dtm/N);
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
        
        function [ fluxS ] = matEvaluateHLLNumFlux1d( obj, nx, fm, fp )
            [ fluxS ] = mxEvaluateHLLNumFlux1d( obj.hmin, obj.gra, nx, fm, fp );
        end
        
        function [ fluxS ] = matEvaluateLFNumFlux1d( obj, mesh, nx, fm, fp )
            
        end
    end
end

