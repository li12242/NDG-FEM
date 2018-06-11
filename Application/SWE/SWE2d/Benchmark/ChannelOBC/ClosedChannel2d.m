classdef ClosedChannel2d < OpenChannel2d
    
    methods
        function obj = ClosedChannel2d( N, M )
            obj = obj@OpenChannel2d( N, M );
        end
        
        function drawGaugePoints( obj )
            waveLen = obj.T * sqrt(obj.gra * obj.H);
            xg = -obj.ChLength:waveLen/6:0;
            drawGaugePoints@OpenChannel2d( obj, xg );
        end
    end
    
    methods( Access = protected, Static )
        function obtype = setOpenBoundaryCondition(  )
            obtype = [ ...
                enumBoundaryCondition.Flather, ...
                enumBoundaryCondition.SlipWall ];
        end
    end
    
    methods( Access = protected )
        function [ option ] = setOption( obj, option )
            ftime = 5000;
            outputIntervalNum = 2000;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = ...
                [mfilename, '.', num2str(obj.meshUnion.cell.N)];
            option('temporalDiscreteType') = enumTemporalDiscrete.RK45;
            option('limiterType') = enumLimiter.None;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('CoriolisType') = enumSWECoriolis.None;
            option('WindType') = enumSWEWind.None;
            option('FrictionType') = enumSWEFriction.None;
        end
        
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                fphys{m}(:, :, 1) = + obj.H;
                fphys{m}(:, :, 4) = - obj.H;
            end
        end
        
        %> Set exact solution for standarding wave
        function [h, hu] = setExactSolution( obj, x, time )
            w = 2 * pi / obj.T;
            c = sqrt( obj.gra * obj.H );
            k = w/c;
            
            temp1 = cos( k .* x ) * cos( w * time );
            temp2 = sin( k .* x ) * sin( w * time );
            eta = obj.delta * obj.H;
            h = 2* eta * temp1 + obj.H;
            u = 2* eta * sqrt( obj.gra/obj.H ) * temp2;
            hu = h .* u;
        end% func
    end
    
end

