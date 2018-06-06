classdef SteadyMount1d < SWETransitionCell1d
    %STEADYMOUNT1D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        hmin = 1e-3
        gra = 9.81
    end
    
    methods
        function obj = SteadyMount1d( N, M )
            obj = obj@SWETransitionCell1d();
            [ mesh ] = makeUniformMesh( N, M );
            obj.initPhysFromOptions( mesh );
        end
    end
    
    methods( Access = protected )
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{1} = zeros(mesh.cell.Np, mesh.K, obj.Nfield);
                
                fphys{1}(:,:,3) = max( 0, 0.25 - 5*(mesh.x - 0.5).^2 );
                fphys{1}(:,:,1) = max( 0, 0.2 - fphys{1}(:,:,3) );
            end
        end
        
        function [ option ] = setOption( obj, option )
            
            ftime = 0.5;
            outputIntervalNum = 50;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('limiterType') = NdgLimiterType.TVB;
            option('limiterParameter') = 1e-5;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
        end
    end
    
end

function [ mesh ] = makeUniformMesh( N, M )
xlim = [0, 1];
bcType = [NdgEdgeType.SlipWall, NdgEdgeType.SlipWall];
[ mesh ] = makeUniformMesh1d( N, xlim, M, bcType );
end

