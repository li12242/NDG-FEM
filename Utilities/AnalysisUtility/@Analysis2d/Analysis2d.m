classdef Analysis2d < AbstractAnalysis
    properties (SetAccess = protected)
        
    end

    methods (Access = public)
        function obj = Analysis2d( phys, xg, yg )
            obj = obj@AbstractAnalysis( phys );
            obj = initGaugePointInfo( obj, xg, yg );
        end
    end

    methods (Access = protected)
        function obj = initGaugePointInfo(obj, xg, yg)
            obj.Ng = numel(xg);
            obj.xg = xg;
            obj.yg = yg;
            [ obj.gaugeMesh, obj.gaugeCell ] = accessGaugePointLocation( obj );
            obj.Vg = assessGaugeInterpMatrix( obj );
        end

        [ gaugeMesh, gaugeCell ] = accessGaugePointLocation( obj );
        [ Vg ] = assessGaugeInterpMatrix( obj );
    end
end