classdef DamBreakWetUniformMesh1d < SWEConventional1d
    
    properties( Constant )
        hmin = 1e-4
        gra = 9.81
        
        x0 = 500
        h0 = 10
        h1 = 2
    end
    
    methods
        function obj = DamBreakWetUniformMesh1d( N, M )
            obj = obj@SWEConventional1d();
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
                ind = (mesh.xc < obj.x0);
                fphys{1}(:,  ind, 1) = obj.h0;
                fphys{1}(:, ~ind, 1) = obj.h1;
            end
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 20;
            outputIntervalNum = 50;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('limiterType') = NdgLimiterType.TVB;
            option('limiterParameter') = 1e-4;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
            option('NumFluxType') = SWENumFluxType.ROE;
        end
        
        function [ fext ] = getExactFunction( obj, time )
            fext = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                fext{m} = zeros(obj.meshUnion(m).cell.Np, obj.meshUnion(m).K, obj.Nfield);
                [ h, u ] = getExtFunc1d(obj, obj.meshUnion(m).x, time);
                hu = h.*u;
                fext{m}(:,:,1) = h;
                fext{m}(:,:,2) = hu;
            end
            
        end
    end
    
end

function [h, u] = getExtFunc1d(obj, x, time)
h = zeros(size(x));
u = zeros(size(x));
temp = (x - obj.x0)/time;
c0 = sqrt(obj.gra .* obj.h0);
% left part
ind = (temp < -c0 );
h(ind) = obj.h0;
u(ind) = 0;
% middle part
ind = (temp > -c0 ) & ( temp < 2*c0 );
h(ind) = (2*c0 - temp(ind) ).^2/obj.gra/obj.gra;
u(ind) = 2/3*(c0 + temp(ind) );
end

function [ mesh ] = makeUniformMesh( N, M )
xlim = [0, 1000];
bcType = [NdgEdgeType.ZeroGrad, NdgEdgeType.ZeroGrad];
[ mesh ] = makeUniformMesh1d( N, xlim, M, bcType );
end


