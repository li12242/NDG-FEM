classdef Malpasset2d < SWEWD2d & SWEPreBlanaced2d
    
    properties( Constant )
        hmin = 1e-1
        n = 0.029.^2
        gra = 9.81
        gmshfile = 'SWE/SWE2d/Benchmark/@Malpasset2d/mesh/malpasset.msh';
    end
    
    
    methods(Access=protected, Hidden)
        function fphys = matInterpolateTopography( obj, fphys ) 
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                
                tpfile = ...
                     'SWE/SWE2d/Benchmark/@Malpasset2d/mesh/bathmetry1.txt';
                fp = fopen(tpfile);
                fgets(fp);
                data = fscanf(fp, '%e %e %e', [3, inf]);
                fclose(fp);
                interp = scatteredInterpolant( ...
                    data(1,:)', data(2,:)', data(3,:)', 'linear');

                fphys{m}(:,:,4) = interp(mesh.x, mesh.y);
            end
        end% func
    end
    
    methods
        function obj = Malpasset2d( N )
            obj = obj@SWEPreBlanaced2d();
            mesh = makeGmshFileUMeshUnion2d( N, obj.gmshfile );
            obj.initPhysFromOptions( mesh );
        end
        
        function assessMexGaugeDepth( obj )
        end
    end
    
    methods(Access=protected)
        function [ option ] = setOption( obj, option )
            ftime = 2000;
            
            outputIntervalNum = 200;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('cfl') = 1/obj.meshUnion(1).cell.N;
            option('temporalDiscreteType') = NdgTemporalIntervalType.DeltaTime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK22;
            option('limiterType') = NdgLimiterType.Vert;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
            option('FrictionType') = FrictionType.Manning;
            option('FrictionCoefficient_n') = obj.n;
        end
        
        function [ fphys ] = matEvaluatePostFunc(obj, fphys)
            [ fphys ] = matEvaluatePostFunc@SWEAbstract2d( obj, fphys );
            obj.matUpdateWetDryState( fphys );
        end% func
        
        function fphys = matEvaluateLimiter( obj, fphys )            
            obj.matUpdateWetDryState( fphys )
            
            fphys = obj.limiter.matLimit( fphys, 2 );
            fphys = obj.limiter.matLimit( fphys, 3 );
            for m = 1:obj.Nmesh % update new elevation
                fphys{m}(:,:,5) = fphys{m}(:,:,1) + fphys{m}(:,:,4);
            end
            fphys = obj.limiter.matLimit( fphys, 5 ); % enforce the elevation
            for m = 1:obj.Nmesh % update new elevation
                fphys{m}(:,:,1) = fphys{m}(:,:,5) - fphys{m}(:,:,4);
            end
            %fphys = matEvaluateLimiter@SWEPreBlanaced2d( obj, fphys );            
        end% func
        
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
            end
                
            fphys = obj.matInterpolateTopography( fphys );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                
                h = 100 - fphys{m}(:, :, 4);
                ind = (mesh.EToR == int8(NdgRegionType.Dry) );
                h(:, ind) = 0;
%                 h(:, ind) = 0 - fphys{m}(:, ind, 4);
%                 h( h < 0 ) = 0;
                fphys{m}(:, :, 1) = h;
            end
        end% func
    end
    
end

