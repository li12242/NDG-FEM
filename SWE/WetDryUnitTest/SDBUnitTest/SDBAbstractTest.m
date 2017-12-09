classdef SDBAbstractTest < SWEAbstractSDB2d
    %SDBABSTRACTTEST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties( Constant )
        hmin = 1e-4
        xlim = [-1, 1]
        ylim = [-1, 1]        
    end
    
    properties( Abstract, Constant )
        zm
        zp
        hm
        hp
    end
    
    methods
        function obj = SDBAbstractTest()
            obj = obj@SWEAbstractSDB2d();
            mesh = makeUniformMesh( obj );
            obj.initPhysFromOptions( mesh );
        end
        
        function evaluateFlux( obj )
            dflux = cell( obj.Nmesh, 1 );
            E = cell( obj.Nmesh, 1 );
            G = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                [ frecM, frecP ] = obj.boundaryHydroReconst( mesh, obj.fphys{m} );
                [ E{m} ] = obj.matEvaluateSurfFlux( mesh, obj.fphys{m} );
                dflux{m} = obj.matEvaluateNumericalFlux ...
                    ( mesh, obj.fphys{m}, obj.fext{m} );
            end
            mesh = obj.meshUnion;
            Nface = obj.meshUnion.cell.Nface;
            Nfp = obj.meshUnion.cell.Nfp(1);
            fL = 2; fR = 4;
            fpL = sub2ind( [Nfp, Nface], 1:Nfp, fL*ones(1, Nfp) );
            fpR = sub2ind( [Nfp, Nface], 1:Nfp, fR*ones(1, Nfp) );
            nxL = mesh.nx(fpL,1); nxR = mesh.nx(fpR,2);
            nyL = mesh.ny(fpL,1); nyR = mesh.ny(fpR,2);
            
            fprintf('-------- numerical flux --------\n');
            fprintf('| Node Id |      L       |       R       |\n');
            for fld = 1:3
                fprintf(' physical field = %d \n', fld);
                for n = 1:mesh.cell.Nfp(1)
                    fprintf('| f at %d |  %f  |   %f   |\n', n, frecM(fpL(n),1, fld), frecP(fpL(n),1, fld));
                    fprintf('| E at %d |  %f  |   %f   |\n', n, E{1}(fpL(n),1, fld), E{1}(fpR(n),2, fld));
                    fprintf('| G at %d |  %f  |   %f   |\n', n, G{1}(fpL(n),1, fld), G{1}(fpR(n),2, fld));
                    numFluxL = - ( dflux{1}(fpL(n),1,fld) ) + nxL(n).*E{1}(fpL(n),1,fld) + nyL(n).*G{1}(fpL(n),1,fld);
                    numFluxR = - ( dflux{1}(fpR(n),2,fld) ) + nxR(n).*E{1}(fpR(n),2,fld) + nyR(n).*G{1}(fpR(n),2,fld);
                    fprintf('| f* at %d |  %f  |  %f  |\n', n, numFluxL, numFluxR);
                    fprintf('| f - f* at %d |  %f  |  %f  |\n', n, ...
                        dflux{1}(fpL(n),1,fld), dflux{1}(fpR(n),2,fld));
                end
            end
        end
    end
    
    methods( Access = protected )
        
        function [ E, G ] = evaluateSurfaceFluxTerm( obj, mesh, fphys )
            [ frecM, ~ ] = obj.boundaryHydroReconst( mesh, fphys );
            
            E = zeros(mesh.cell.TNfp, mesh.K, 4);
            G = zeros(mesh.cell.TNfp, mesh.K, 4);

            [ um, vm ] = obj.evaluateVelocity(...
                frecM(:,:,1), frecM(:,:,2), frecM(:,:,3) );

            E(:,:,1) = frecM(:,:,2);
            G(:,:,1) = frecM(:,:,3);
            E(:,:,2) = 0.5*obj.gra*frecM(:,:,1).^2 + frecM(:,:,2) .* um;
            G(:,:,2) = frecM(:,:,2) .* vm;
            E(:,:,3) = frecM(:,:,3) .* um;
            G(:,:,3) = 0.5*obj.gra*frecM(:,:,1).^2 + frecM(:,:,3) .* vm;
        end
        
        function [ frecM, frecP ] = boundaryHydroReconst( obj, mesh, fphys )

            eidM = mesh.eidM + mesh.cell.TNfp * 3;
            eidP = mesh.eidP + mesh.cell.TNfp * 3;

            bm = fphys(eidM); bp = fphys(eidP);
            zs = max( bm, bp );
            etam = bm + fphys( mesh.eidM );
            etap = bp + fphys( mesh.eidP );
            zs = min( zs, etam );
            hm = etam - zs;
            hp = max( 0, etap - zs ) - max( 0, bp - zs );

            [ um, vm ] = obj.evaluateVelocity( fphys( mesh.eidM ),...
               fphys( mesh.eidM + mesh.cell.TNfp ), ...
               fphys( mesh.eidM + mesh.cell.TNfp*2 ) );
            [ up, vp ] = obj.evaluateVelocity( fphys( mesh.eidP ),...
               fphys( mesh.eidP + mesh.cell.TNfp ), ...
               fphys( mesh.eidP + mesh.cell.TNfp*2 ) );

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
            ind = ( h>obj.hmin );
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
        end
    end
end

function [ mesh ] = makeUniformMesh( test )
bctype = [...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad];

mesh = makeUniformQuadMesh(1, test.xlim, test.ylim, 2, 1, bctype);
end% func
