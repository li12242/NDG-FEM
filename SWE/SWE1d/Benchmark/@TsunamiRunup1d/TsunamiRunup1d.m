classdef TsunamiRunup1d < SWEWD1d
    %TSUNAMIRUNUP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        hmin = 5e-2;
        gra = 9.81;
        M = 1e-2;
    end
    
    methods
        function obj = TsunamiRunup1d( N, M )
            obj = obj@SWEWD1d();
            mesh = makeUniformMesh( N, M );
            obj.initPhysFromOptions( mesh );
            obj.fext = obj.fphys;
        end
    end
    
    methods( Access = protected )
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                % initial condition
                fphys{m} = zeros(mesh.cell.Np, mesh.K, obj.Nfield);
                % bottom topography
                bot = 5000 - 0.1*mesh.x;
                [path, ~, ~] = fileparts(mfilename('fullpath'));
                data = load([path, '/initial_condition.mat']);
                Interp = griddedInterpolant(data.x, data.eta+5000, 'nearest');
                eta = Interp( mesh.x(:) );
                eta = reshape( eta, mesh.cell.Np, mesh.K );
                h = max( 0, eta - bot) ;
                fphys{m}(:,:,1) = h;
                fphys{m}(:,:,3) = bot;
            end
        end
        
        function [ option ] = setOption( obj, option )            
            ftime = 360;
            outputIntervalNum = 80;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('cfl') = 1/obj.meshUnion(1).cell.N;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('limiterType') = NdgLimiterType.TVB;
            option('limiterParameter') = 1e-10;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
        end
    end
end

function [ mesh ] = makeUniformMesh( N, M )
xlim = [-500, 50000];
bcType = [NdgEdgeType.Inner, NdgEdgeType.ClampedVel];
[ mesh ] = makeUniformMesh1d( N, xlim, M, bcType );
end
