classdef TsunamiRunup2d < SWEPreBlanaced2d

    properties(Constant)
        gra = 9.81
        hmin = 1e-3
        n = 0.01.^2
    end

    properties
        boundaryInfo
    end

    methods (Access = public)
        function obj = TsuamiRunup2d( N )
            [ path, ~, ~ ] = fileparts( mfilename('fullpath') );
            meshfile = [ path, '/mesh/quad.msh' ];
            mesh = makeGmshFileUMeshUnion2d( N, meshfile );
            obj.initPhysFromOptions( mesh );

            % initialize boundary 
            obj.boundaryInfo = obj.initBoundaryInformation( );
        end
        
        %> Compared numerical water elevation with measured data
        checkGaugePoints( obj );
    end

    methods (Access = protected)
        function matUpdateExternalField( obj, time, fphys )
            % find previous input boundary time step
            [ idP ] = find( obj.boundaryInfo.time < time, 1, 'last' );
            % find next input boundary time step
            [ idN ] = find( obj.boundaryInfo.time > time, 1, 'first' );
            if isempty( idP )
                idP = 1;
                delta_time_P = 0;
            else
                delta_time_P = abs( time - obj.boundaryInfo.time(idP) );
            end
            
            if isempty( idN )
                idN = numel( obj.boundaryInfo.time );
                delta_time_N = 0;
            else
                delta_time_N = abs( time - obj.boundaryInfo.time(idN) );
            end
            
            alpha_P = delta_time_P / ( delta_time_P + delta_time_N );
            alpha_N = delta_time_N / ( delta_time_P + delta_time_N );
            eta = obj.boundaryInfo.extElevation(idP) * alpha_P + ...
                obj.boundaryInfo.extElevation(idN) * alpha_N;
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                ind = obj.boundaryInfo.nodeId{m};
                obj.fext{m}( ind ) = eta - fphys{m}( ind + 3 * mesh.cell.Np * mesh.K );
            end
        end

        function [ option ] = setOption( obj, option )
            ftime = 22.5;
            
            outputIntervalNum = 200;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK22;
            option('limiterType') = NdgLimiterType.Vert;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
            option('FrictionType') = SWEFrictionType.Manning;
            option('FrictionCoefficient_n') = obj.n;
        end

        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                fphys{m}(:, :, 4) = obj.interpBottomLevel( mesh );
                h = - fphys{m}(:, :, 4);
                h(h<0) = 0; 
                fphys{m}(:, :, 1) = h;
            end
        end
    end

    methods( Access = private )
        function boundaryInfo = initBoundaryInformation( obj )
            boundaryInfo = struct('nodeId', {}, 'time', {}, 'extElevation', {});
            
            nodeId = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                ind = ( mesh.eidtype == int8( NdgEdgeType.ClampedDepth ) );
                nodeId{m} = mesh.eidM(ind);
            end
            [ path, ~, ~ ] = fileparts( mfilename('fullpath') ); 
            obcfile = [ path, '/mesh/Benchmark_2_input.txt'];
            fp = fopen(obcfile);
            fgetl(fp); % pass the rest of first line
            data = fscanf(fp, '%f %f', [2, inf]);
            fclose(fp);
            
            boundaryInfo(1).nodeId = nodeId;
            boundaryInfo(1).time = data(1, :);
            boundaryInfo(1).extElevation = data(2, :);
        end
    end

    methods(Access = private, Static)
        function bot = interpBottomLevel( mesh )
            [ path, ~, ~ ] = fileparts( mfilename('fullpath') ); 
            topography_file = [ path, ...
                '/mesh/Benchmark_2_Bathymetry.txt'];
            fp = fopen(topography_file);
            fgetl(fp);
            data = fscanf(fp, '%e %e %e', [3, inf]);
            fclose(fp);
            interp = scatteredInterpolant( ...
                data(1,:)', data(2,:)', -data(3,:)','linear');
            bot = interp( mesh.x, mesh.y );
        end
    end
end
