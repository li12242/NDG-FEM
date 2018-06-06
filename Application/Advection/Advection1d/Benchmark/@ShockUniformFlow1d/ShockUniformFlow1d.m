classdef ShockUniformFlow1d < AdvAbstractConstFlow1d
    
    properties
        u0 = 0.5
        x0 = 0.2
        N
    end
    
    methods
        function obj = ShockUniformFlow1d( N, M )
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
                fphys{m} = obj.getExtFunc( obj.meshUnion(m), 0 );
            end
        end
        
        function option = setOption( obj, option )
            outputIntervalNum = 50;
            option('startTime') = 0.0;
            option('finalTime') = 1.2;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = 2.4/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK22;
            
            dt = nan;
            for m = 1:obj.Nmesh
                dt = min( dt, ...
                    min( obj.meshUnion(m).LAV )./obj.u0./(2*obj.N + 1) );
            end
            option('timeInterval') = dt;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
            option('limiterType') = NdgLimiterType.TVB;
            option('limiterParameter') = 5e2;
        end
        
        function [ fext ] = getExtFunc( obj, mesh, time )
            xt = obj.x0 + time * obj.u0;
            ind = abs( mesh.x - xt ) < 0.1;
            fext = zeros( size( mesh.x ) );
            fext(ind) = 1;
        end
        
    end
    
end

function [ mesh ] = makeUniformMesh( N, M )
xlim = [0, 1];
bcType = [NdgEdgeType.Clamped, NdgEdgeType.Clamped];
[ mesh ] = makeUniformMesh1d( N, xlim, M, bcType );
end

