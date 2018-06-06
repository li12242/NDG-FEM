classdef AbstractAnalysis < handle
    %ABSTRACTPOSTPROCESS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        phys
    end
    
    properties ( SetAccess = protected )
        %> num of gauge points
        Ng
        %> gauge point position
        xg, yg, zg
        %> mesh id of gauge point location
        gaugeMesh
        %> cell id of gauge point location
        gaugeCell
        %> interp matrix
        Vg
    end

    methods (Abstract, Access = protected)
        [ gaugeMesh, gaugeCell ] = accessGaugePointLocation( obj );
        [ Vg ] = assessGaugeInterpMatrix( obj );
    end
    
    methods
        function obj = AbstractAnalysis( phys )
            obj.phys = phys;
        end

        InterpGaugeResult( obj, fieldId )
    end
    
end

