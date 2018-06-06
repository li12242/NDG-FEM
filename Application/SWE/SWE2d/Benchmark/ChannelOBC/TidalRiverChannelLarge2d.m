classdef TidalRiverChannelLarge2d < TidalRiverChannel2d
    
    properties( Constant )
        %> left position of channel
        ChLeft = - 28e3 * 3; % domain length = wave length * lambda
        %> right position of channel
        ChRight = 28e3 * 6;
    end
    
    methods( Access = public )
        function obj = TidalRiverChannelLarge2d( N, M )
            obj = obj@TidalRiverChannel2d( N, M );
        end
    end
    
    methods( Access = protected )
        function [ option ] = setOption( obj, option )
            ftime = 1.3 * ( obj.ChRight - obj.ChLeft ) / sqrt( obj.gra * obj.H );
            outputIntervalNum = 1000;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = ...
                [mfilename, '.', num2str(obj.meshUnion.cell.N)];
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK22;
            option('limiterType') = NdgLimiterType.None;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.GaussQuadrature;
            option('CoriolisType') = SWECoriolisType.None;
            option('WindType') = SWEWindType.None;
            option('FrictionType') = SWEFrictionType.None;
        end
    end
    
    methods( Access = protected, Static  )
        %> set open boundary condition
        function obtype = setOpenBoundaryCondition(  )
            obtype = [NdgEdgeType.Clamped, NdgEdgeType.ClampedVel];
        end
    end
    
end

