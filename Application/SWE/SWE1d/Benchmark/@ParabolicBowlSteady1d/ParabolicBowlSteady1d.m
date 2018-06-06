classdef ParabolicBowlSteady1d < ParabolicBowl1d
    
    properties
        eta0 = 0.1
    end
    
    methods
        function obj = ParabolicBowlSteady1d( N, M )
            obj = obj@ParabolicBowl1d(N, M);
        end
    end
    
    methods( Access = protected )
        function fphys = setInitialField( obj )
            fphys = obj.evaluateExactFunc(0);
        end
        
        function fphys = evaluateExactFunc( obj, time )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                
                fphys{m} = zeros(mesh.cell.Np, mesh.K, obj.Nfield);
                % Init Condition
                fphys{m}(:,:,3) = obj.h0.*(mesh.x.^2./obj.a^2 - 1);
                z = obj.eta0;
                fphys{m}(:,:,1) = max( 0, z - fphys{m}(:,:,3));
            end
        end
        
        function [ fphys ] = matEvaluatePostFunc(obj, fphys)
            fphys = matEvaluatePostFunc@ParabolicBowl1d(obj, fphys);
            fphys{1}(:, :, 2) = fphys{1}(:,:,2) * 0.98;
        end
        
        function [ option ] = setOption( obj, option )
            
            %T = 2*pi*obj.a/sqrt(2*obj.gra*obj.h0);
            %ftime = 2*T;
            ftime = 2e3;
            outputIntervalNum = 1000;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('cfl') = 1/obj.meshUnion(1).cell.N;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('limiterType') = NdgLimiterType.TVB;
            option('limiterParameter') = 1;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
        end
    end
    
end

