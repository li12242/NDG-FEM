classdef CSBAbstractTest < SWEPreBlanaced2d
    
    methods
        function checkWellBalanced( obj, m, k )
            mesh = obj.meshUnion(m);
            Np = mesh.cell.Np;
            % variable
            h = obj.fphys{m}(:, k, 1);
            hu = obj.fphys{m}(:, k, 2);
            hv = obj.fphys{m}(:, k, 3);
            z = obj.fphys{m}(:, k, 4);
            fprintf('\nConservative Variable\n');
            for n = 1:Np
                fprintf('n = %d, h = %f, hu = %f, hv = %f, z = %f, eta = %f\n', ...
                    n, h(n), hu(n), hv(n), z(n), h(n) + z(n));
            end
            
            % volume integral
            fprintf('\nCheck Volume Flux Term\n');
            obj.advectionSolver.evaluateAdvectionRHS( obj.fphys );
            [ E, G ] = obj.matEvaluateFlux( mesh, obj.fphys{m} );
            eflux = E(:, k, :);
            gflux = G(:, k, :);
            fprintf('Flux Term\n');
            for n = 1:Np
                fprintf('n = %d, E(h) = %f, E(hu) = %f, E(hv) = %f\n', ...
                    n, eflux(n, :, 1), eflux(n, :, 2), eflux(n, :, 3));
            end
            
            for n = 1:Np
                fprintf('n = %d, G(h) = %f, G(hu) = %f, G(hv) = %f\n', ...
                    n, gflux(n, :, 1), gflux(n, :, 2), gflux(n, :, 3));
            end
            
            fprintf('\nExact Volume Flux Term\n');
            for n = 1:Np
                [ Eh, Ehu, Ehv, Gh, Ghu, Ghv ] = obj.evaluateVolumeFluxTerm...
                    ( h(n), hu(n), hv(n), z(n) );
                
                fprintf('n = %d, E(h) = %f, E(hu) = %f, E(hv) = %f\n', ...
                    n, Eh, Ehu, Ehv );
                fprintf('n = %d, G(h) = %f, G(hu) = %f, G(hv) = %f\n', ...
                    n, Gh, Ghu, Ghv );
            end
%             fprintf('|  Field  |  Node Id  |  E  |  G  |  F  |\n');
%             for fld = 1:3
%                 for n = 1:cell.Np
%                     fprintf('|  fld = %d  |  %d  |  E = %e  |  G = %e  |  F = %e  |\n', fld, n, ...
%                         E(n, k, fld), G(n, k, fld), obj.frhs{m}(n, k, fld)  );
%                 end
%             end
            % volume integral
%             obj.frhs{m} = zeros( obj.meshUnion.cell.Np, obj.meshUnion.K );
%             obj.matEvaluateSourceTerm( obj.fphys );
%             fprintf('-------- topography source --------\n');
%             for fld = 1:3
%                 for n = 1:cell.Np
%                     fprintf('|  fld = %d  | %d | S = %e  |\n', fld, n, obj.frhs{m}(n, k, fld) );
%                 end
%             end
%             % rhs
%             obj.matEvaluateRHS( obj.fphys );
%             fprintf('-------- RHS --------\n');
%             for fld = 1:3
%                 for n = 1:cell.Np
%                     fprintf('|  fld = %d  | %d | rhs = %e  |\n', fld, n, obj.frhs{m}(n, k, fld) );
%                 end
%             end
        end% func
    end
    
    methods( Access = protected )
        function [ Eh, Ehu, Ehv, Gh, Ghu, Ghv ] = evaluateVolumeFluxTerm( obj, h, hu, hv, z )
            Eh = hu;
            Gh = hv;
            Ehu = ( hu.^2 / h + 0.5 * obj.gra * (h ^ 2 - z ^ 2) );
            Ehv = hu * hv / h;
            Ghu = hu * hv / h;
            Ghv = ( hv.^2 / h + 0.5 * obj.gra * (h ^ 2 - z ^ 2) );
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
            option('limiterType') = NdgLimiterType.Vert;
        end
    end
end


