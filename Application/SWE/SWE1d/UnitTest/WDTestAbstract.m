classdef WDTestAbstract < SWEWDMesh1d
    %WDTEST1 Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = WDTestAbstract()
        end
        
        function TestWetDryReconstruct( obj )
            obj.fphys = obj.matEvaluatePostFunc( obj.fphys );
        end
        
        function checkVolumTerm( obj, m, k1 )
            mesh = obj.meshUnion(m);
            [ E ] = obj.matEvaluateFlux( mesh, obj.fphys{m} );
            cell = mesh.cell;
            fprintf('-------- volume term check routine --------\n');
            fprintf('m1 = %d, f1 = %d\n', m, k1);
            fprintf('Region type = \n');
            NdgRegionType( mesh.EToR(k1) )
            for n = 1:cell.Np
                fprintf('n = %d, h = %f, qx = %f, z = %f\n', n, ...
                    obj.fphys{m}(n, k1, 1), ...
                    obj.fphys{m}(n, k1, 2), ...
                    obj.fphys{m}(n, k1, 3) );
            end
            for n = 1:cell.Np
                fprintf('n = %d, E(h) = %f, E(qx) = %f\n', n, ...
                    E(n, k1, 1), ...
                    E(n, k1, 2) );
            end
            
            obj.advectionSolver.evaluateAdvectionRHS( obj.fphys );
            for n = 1:cell.Np
                fprintf('n = %d, rhs(h) = %f, rhs(qx) = %f\n', n, ...
                    obj.frhs{m}(n, k1, 1), ...
                    obj.frhs{m}(n, k1, 2) );
            end
        end
        
        function checkFluxTerm( obj, m, k1, f1 )
            obj.matUpdateExternalField( 0, obj.fphys );
            mesh = obj.meshUnion(m);
            nx = obj.meshUnion.nx;
            [ fM, fP ] = obj.matEvaluateSurfaceValue( mesh, obj.fphys{m}, obj.fext{m} );
            [ flux  ] = obj.matEvaluateSurfFlux( mesh, nx, fM );
            [ fluxS ] = obj.matEvaluateSurfNumFlux( mesh, nx, fM, fP );
            
            fprintf('-------- numerical flux check routine --------\n');
            fprintf('k1 = %d, f1 = %d\n', k1, f1);
            fprintf('Adjacent cell:\n')
            k2 = mesh.EToE( f1, k1 );
            f2 = mesh.EToF( f1, k1 );
            fprintf('k1 = %d, f1 = %d\n', k2, f2);
            fid = sum( mesh.cell.Nfp(1:f1) - mesh.cell.Nfp(f1) + 1 );
            fprintf('-------- Riemann State at Boundary --------\n');
            fprintf('edge type = \n');
            NdgEdgeType( mesh.eidtype( fid, k1 ) )
            fprintf('|    | \tU1\t | \tU2\t |\n');
            fprintf('| h  | \t%f\t | \t%f\t |\n', fM(fid, k1, 1), fP(fid, k1, 1));
            fprintf('| qx | \t%f\t | \t%f\t |\n', fM(fid, k1, 2), fP(fid, k1, 2));
            fprintf('| F* | \t%f\t | \t%f\t |\n', fluxS(fid, k1, 2), fluxS(fid, k1, 2));
            fprintf('| F- | \t%f\t | \t%f\t |\n', flux(fid, k1, 2), flux(fid, k1, 2));
        end
    end
    
    methods( Access = protected )
        function [ option ] = setOption( obj, option )
            ftime = 20;
            outputIntervalNum = 50;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('limiterType') = NdgLimiterType.None;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
        end
    end
    
end



