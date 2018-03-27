classdef OpenChannel2d < SWEWD2d
    
    properties(Constant)
        hmin = 1e-4; % water depth threshold
        gra = 9.81; % gravity acceleration
        %> channel depth
        H = 40; 
        %> amplitude
        eta = 0.2;
        %> wave period 6 min
        T = 360;

        %> domain length
        ChLength = 2e3;
        %> domain width
        ChWidth = 500;
        %> mesh rotation angle
        theta = 0;
        %> mesh center 
        xc = 0;
        %> mesh center 
        yc = 0;
    end
    
    methods
        function obj = OpenChannel2d( N, type )
            obj = obj@SWEWD2d();
            mesh = obj.makeUniformMesh( N, type );
            obj.initPhysFromOptions( mesh );
        end
    end

    methods( Access = protected, Hidden )
        %> create uniform computation grid
        makeUniformMesh( obj, N, type )
    end

    methods( Access = protected )
        function [ option ] = setOption( obj, option )
            ftime = 7200;
            outputIntervalNum = 2000;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('temporalDiscreteType') = NdgTemporalIntervalType.DeltaTime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('limiterType') = NdgLimiterType.None;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
            option('CoriolisType') = SWECoriolisType.None;
            option('WindType') = SWEWindType.None;
            option('FrictionType') = SWEFrictionType.None;
        end

        function fphys = setInitialField( obj )
            w = 2 * pi / obj.T;
            c = sqrt( obj.gra * obj.H );
            k = w/c;

            fphys = cell( obj.Nmesh );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                h = obj.eta * cos( k * mesh.x) + obj.H;
                u = obj.eta * sqrt( obj.gra /obj.H ) * cos( k * obj.mesh.x );
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                fphys{m}(:, :, 1) = h;
                fphys{m}(:, :, 2) = h .* u;
            end
        end% func
    end% methods
    
end

