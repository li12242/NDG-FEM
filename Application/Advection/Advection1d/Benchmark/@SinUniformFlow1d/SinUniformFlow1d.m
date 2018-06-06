classdef SinUniformFlow1d < AdvAbstractConstFlow1d
    
    properties
        u0 = 0.5
        N
    end
    
    methods
        function obj = SinUniformFlow1d( N, M )
            obj = obj@AdvAbstractConstFlow1d();
            obj.N = N;
            [ mesh ] = makeUniformMesh( N, M );
            obj.initPhysFromOptions( mesh );
        end        
    end
    
    methods( Access = protected )
        function [ fphys ] = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                fphys{m} = sin( obj.meshUnion(m).x );
            end
        end
        
        function option = setOption( obj, option )
            outputIntervalNum = 50;
            option('startTime') = 0.0;
            option('finalTime') = 24;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = 2.4/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            
            dt = nan;
            for m = 1:obj.Nmesh
                dt = min( dt, ...
                    min( obj.meshUnion(m).LAV )./obj.u0./(2*obj.N + 1) );
            end
            option('timeInterval') = dt;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
            option('limiterType') = NdgLimiterType.None;
        end
        
        function [ fext ] = getExtFunc( obj, mesh, time )
            fext = sin( mesh.x - time * obj.u0 );
        end
        
    end
    
end

function [ mesh ] = makeUniformMesh( N, M )
xlim = [0, 2*pi];
bcType = [NdgEdgeType.Clamped, NdgEdgeType.Clamped];
[ mesh ] = makeUniformMesh1d( N, xlim, M, bcType );
end
