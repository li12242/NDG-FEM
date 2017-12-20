classdef CSBAbstractTest < SWEPreBlanaced2d
    
    properties( Constant )
%         hmin = 1e-4
%         gra = 9.81
%         xlim = [-1, 1];
%         ylim = [-1, 1];
%         h0 = 2;
    end
    
    properties( Abstract, Constant )
%         bx
%         by
    end
    
    methods
        function obj = CSBAbstractTest()
        end
        
        function checkWellBalanced( obj, m, k )
            mesh = obj.meshUnion(m);
            % volume integral
            obj.advectionSolver.evaluateAdvectionRHS( obj.fphys );
            [ E, G ] = obj.matEvaluateFlux( mesh, obj.fphys{m} );
            cell = mesh.cell;
            fprintf('-------- volume flux --------\n');
            fprintf('|  Field  |  Node Id  |  E  |  G  |  F  |\n');
            for fld = 1:3
                for n = 1:cell.Np
                    fprintf('|  fld = %d  |  %d  |  E = %e  |  G = %e  |  F = %e  |\n', fld, n, ...
                        E(n, k, fld), G(n, k, fld), obj.frhs{m}(n, k, fld)  );
                end
            end
            % volume integral
            obj.frhs{m} = zeros( obj.meshUnion.cell.Np, obj.meshUnion.K );
            obj.matEvaluateSourceTerm( obj.fphys );
            fprintf('-------- topography source --------\n');
            for fld = 1:3
                for n = 1:cell.Np
                    fprintf('|  fld = %d  | %d | S = %e  |\n', fld, n, obj.frhs{m}(n, k, fld) );
                end
            end
            % rhs
            obj.matEvaluateRHS( obj.fphys );
            fprintf('-------- RHS --------\n');
            for fld = 1:3
                for n = 1:cell.Np
                    fprintf('|  fld = %d  | %d | rhs = %e  |\n', fld, n, obj.frhs{m}(n, k, fld) );
                end
            end
        end
    end
    
    methods( Access = protected )
        
        function fphys = setInitialField( obj )
            fphys = cell(1, 1);
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{1} = zeros( obj.meshUnion.cell.Np, obj.meshUnion.K, obj.Nfield );
                fphys{1}(:,:,4) = obj.bx * mesh.x + obj.by * mesh.y;
                fphys{1}(:,:,1) = obj.h0 - fphys{1}(:,:,4);
            end
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 0.2;
            outputIntervalNum = 2;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('temporalDiscreteType') = NdgTemporalIntervalType.DeltaTime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('WellBlancedType') = true;
        end
    end
end


