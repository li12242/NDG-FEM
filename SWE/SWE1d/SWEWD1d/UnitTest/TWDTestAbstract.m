classdef TWDTestAbstract < SWETransitionCell1d
    %WDTEST1 Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = TWDTestAbstract()
        end
    end
    
    methods( Access = protected )
        function [ option ] = setOption( obj, option )
            ftime = 1;
            outputIntervalNum = 50;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('limiterType') = NdgLimiterType.None;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
        end
    end
    
end
