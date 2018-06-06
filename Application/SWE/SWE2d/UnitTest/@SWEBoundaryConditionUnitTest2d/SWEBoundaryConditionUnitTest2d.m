classdef SWEBoundaryConditionUnitTest2d < SWEConventional2d
    
    properties( Constant )
        gra = 9.81
        hmin = 1e-4
    end
    
    methods
        function obj = SWEBoundaryConditionUnitTest2d( N, bctype )
            mesh = makeUniformQuadMesh(N, [0, 10], [0, 10], 2, 2, bctype);
            obj.initPhysFromOptions( mesh );
        end
        
        function checkBoundaryValue( obj, k )
            obj.matUpdateExternalField( 0, obj.fphys );
            [ bctype ] = obj.meshUnion.EToB(:, k);
            for f = 1:numel(bctype)
                if ( bctype(f) ~= NdgEdgeType.Inner )
                    obj.printBoundaryCondition( k, f, ...
                        obj.fphys{1}(:,:, 1), ...
                        obj.fphys{1}(:,:, 2), ...
                        obj.fphys{1}(:,:, 3), ...
                        obj.fphys{1}(:,:, 4), ...
                        obj.fext{1}(:, :, 1), ...
                        obj.fext{1}(:, :, 2), ...
                        obj.fext{1}(:, :, 3), ...
                        obj.fext{1}(:, :, 4) )
                end
            end
        end
    end
    
    methods( Access = private )
        function printBoundaryCondition( obj, k, f, h, hu, hv, z, he, hue, hve, ze )
            Nfp = obj.meshUnion.cell.Nfp(f);
            FFP = sum( obj.meshUnion.cell.Nfp(1:f) );
            eidF = FFP:-1:(FFP - obj.meshUnion.cell.Nfp(f)+1 );
            eidM = obj.meshUnion.eidM( eidF, k );
            fprintf('Inner physical field at boundary: %d\n', f);
            for n = 1:Nfp
                ind = eidM(n);
                fprintf('n=%d, h = %f, hu = %f, hv = %f, z = %f, eta = %f\n', ...
                    n, h(ind), hu(ind), hv(ind), z(ind), h(ind) + z(ind) );
            end
            
            fprintf('External field at boundary: %d\n', f);
            for n = 1:Nfp
                ind = eidM(n);
                fprintf('n=%d, h = %f, hu = %f, hv = %f, z = %f, eta = %f\n', ...
                    n, he(ind), hue(ind), hve(ind), ze(ind), he(ind) + ze(ind) );
            end
            
            [fM, fP] = obj.matEvaluateSurfaceValue( obj.meshUnion, obj.fphys{1}, obj.fext{1} );
            Ntol = obj.meshUnion.K * obj.meshUnion.cell.TNfp;
            fprintf('Inner field at boundary: %d\n', f);
            for n = 1:Nfp
                ind = eidF(n) + (k-1) * obj.meshUnion.cell.TNfp;
                gind = eidM(n);
                fprintf('n=%d, h = %f, hu = %f, hv = %f, z = %f, eta = %f\n', ...
                    n, ...
                    fM(ind), ...
                    fM(ind + Ntol), ...
                    fM(ind + 2*Ntol), ...
                    z(gind), ...
                    fM(ind) + z(gind) );
            end
            
            fprintf('Extern field at boundary: %d\n', f);
            for n = 1:Nfp
                ind = eidF(n) + (k-1) * obj.meshUnion.cell.TNfp;
                gind = eidM(n);
                fprintf('n=%d, h = %f, hu = %f, hv = %f, z = %f, eta = %f\n', ...
                    n, ...
                    fP(ind), ...
                    fP(ind + Ntol), ...
                    fP(ind + 2*Ntol), ...
                    z(gind), ...
                    fP(ind) + z(gind) );
            end
        end
    end
    
    methods( Access = protected )
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                fphys{m}(:,:,1) = 2;
                fphys{m}(:,:,2) = 0.5;
                fphys{m}(:,:,4) = -1;
            end
        end
        
        function matUpdateExternalField( obj, time, fphys )
            for m = 1:obj.Nmesh
                obj.fext{m}(:, :, 1) = 5;
                obj.fext{m}(:, :, 2) = 0.2;
                obj.fext{m}(:, :, 3) = 0.3;
            end
        end
        
        function option = setOption( obj, option )
            ftime = 2;
            outputIntervalNum = 2;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('temporalDiscreteType') = NdgTemporalIntervalType.DeltaTime;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = [mfilename, '.', num2str(obj.meshUnion.cell.N)];
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('limiterType') = NdgLimiterType.None;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.GaussQuadrature;
            option('CoriolisType') = SWECoriolisType.None;
            option('WindType') = SWEWindType.None;
            option('FrictionType') = SWEFrictionType.None;
        end
    end
    
end

