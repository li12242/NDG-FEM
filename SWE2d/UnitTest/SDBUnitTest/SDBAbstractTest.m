classdef SDBAbstractTest < SWEPreBlanaced2d
    
    methods
        function obj = SDBAbstractTest
        end
        
        function checkFluxTerm( obj, m, k1, f1 )
            
            mesh = obj.meshUnion(m);
            nx = obj.meshUnion.nx;
            ny = obj.meshUnion.ny;
            [ fM, fP ] = obj.matEvaluateSurfaceValue( mesh, obj.fphys{m}, obj.fext{m} );
            [ frecM, frecP ] = obj.boundaryHydroReconst( mesh, obj.fphys{m} );
            [ flux  ] = obj.matEvaluateSurfFlux( mesh, nx, ny, fM );
            [ fluxS ] = obj.matEvaluateSurfNumFlux( mesh, nx, ny, fM, fP );

            K = [k1, 2];
            fL = f1;
            K(2) = mesh.EToE( fL, K(1) );
            fR = mesh.EToF( fL, K(1) );
            efp = sum( mesh.cell.Nfp(1:fL) );
            sfp = efp - mesh.cell.Nfp(fL) + 1;
            fpm = sfp:efp;

            efp = sum( mesh.cell.Nfp(1:fR) );
            sfp = efp - mesh.cell.Nfp(fR) + 1;
            fpp = sfp:efp;

            fprintf('-------- numerical flux --------\n');
            fprintf('| Node Id |      L       |       R       |\n');
            fprintf('| z- and z+ at %d |  %f  |   %f   |\n', K(1), frecM( fpm(1), K(1), 4), frecP( fpm(1), K(1), 4) );
            fprintf('| z+ and z- at %d |  %f  |   %f   |\n', K(2), frecM( fpp(1), K(2), 4), frecP( fpp(1), K(2), 4) );

            for fld = 1:3
                fprintf(' physical field = %d \n', fld);
                fprintf('| U- and U+ at %d |  %f  |   %f   |\n', K(1), fM(fpm(1), K(1), fld), fP(fpm(1), K(1), fld) );
                fprintf('| Ue- and Ue+ at %d |  %f  |   %f   |\n', K(1), frecM(fpm(1), K(1), fld), frecP(fpm(1), K(1), fld) );
                fprintf('| U+ and U- at %d |  %f  |   %f   |\n', K(2), fM(fpp(1), K(2), fld), fP(fpp(1), K(2), fld) );
                fprintf('| Ue+ and Ue- at %d |  %f  |   %f   |\n', K(2), frecM(fpp(1), K(2), fld), frecP(fpp(1), K(2), fld) );
                fprintf('| F- and F+ at %d |  %f  |   %f   |\n', K(1), flux(fpm(1), K(1), fld), flux(fpp(1), K(2), fld) );
                fprintf('| F+ and F- at %d |  %f  |   %f   |\n', K(2), flux(fpp(1), K(2), fld), flux(fpm(1), K(1), fld) );
                fprintf('| Flux* at %d |  %f  |  %f  |\n', K(1), fluxS(fpm(1), K(1), fld), fluxS(fpp(1), K(2), fld) );
            end
            
        end
    end
    
    methods( Access = protected )
        
        function [ frecM, frecP ] = boundaryHydroReconst( obj, mesh, fphys )
            
            eidM = mesh.eidM + mesh.cell.Np*mesh.K * 3;
            eidP = mesh.eidP + mesh.cell.Np*mesh.K * 3;
            
            bm = fphys( eidM );
            bp = fphys( eidP );
            zs = max( bm, bp );
            etam = bm + fphys( mesh.eidM );
            etap = bp + fphys( mesh.eidP );
            zs = min( zs, etam );
            hm = etam - zs;
            hp = max( 0, etap - zs ) - max( 0, bp - zs );
            
            [ um, vm ] = obj.evaluateVelocity( fphys( mesh.eidM ),...
                fphys( mesh.eidM + mesh.cell.Np*mesh.K ), ...
                fphys( mesh.eidM + mesh.cell.Np*mesh.K*2 ) );
            [ up, vp ] = obj.evaluateVelocity( fphys( mesh.eidP ),...
                fphys( mesh.eidP + mesh.cell.Np*mesh.K ), ...
                fphys( mesh.eidP + mesh.cell.Np*mesh.K*2 ) );
            
            hum = hm .* um; hup = hp .* up;
            hvm = hm .* vm; hvp = hp .* vp;
            
            frecM = zeros(mesh.cell.TNfp, mesh.K, 4);
            frecP = zeros(mesh.cell.TNfp, mesh.K, 4);
            frecM(:,:,1) = hm; frecM(:,:,2) = hum;
            frecP(:,:,1) = hp; frecP(:,:,2) = hup;
            frecM(:,:,4) = zs; frecM(:,:,3) = hvm;
            frecP(:,:,4) = zs; frecP(:,:,3) = hvp;
        end
        
        function [ u, v ] = evaluateVelocity( obj, h, hu, hv )
            u = zeros( size(h) );
            v = zeros( size(h) );
            ind = ( h > obj.hmin );
            u(ind) = hu(ind) ./ h(ind);
            v(ind) = hv(ind) ./ h(ind);
        end
        
        function fphys = setInitialField( obj )
            fphys = cell(1, 1);
            fphys{1} = zeros( obj.meshUnion.cell.Np, obj.meshUnion.K, obj.Nfield );
            fphys{1}(:,1,1) = obj.hm;
            fphys{1}(:,2,1) = obj.hp;
            fphys{1}(:,1,4) = obj.zm;
            fphys{1}(:,2,4) = obj.zp;
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
    
    methods( Static )
        function [ mesh ] = makeUniformMesh( xlim, ylim )
            bctype = [...
                NdgEdgeType.ZeroGrad, NdgEdgeType.ZeroGrad, ...
                NdgEdgeType.ZeroGrad, NdgEdgeType.ZeroGrad];
            
            mesh = makeUniformQuadMesh(1, xlim, ylim, 2, 1, bctype);
        end% func
    end
end


