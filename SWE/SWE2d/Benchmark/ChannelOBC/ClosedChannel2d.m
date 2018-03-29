classdef ClosedChannel2d < OpenChannel2d
    
    methods
        function obj = ClosedChannel2d( N, M, type )
            obj = obj@OpenChannel2d( N, M, type );
        end
        
        function drawGaugePoints( obj )
            waveLen = obj.T * sqrt(obj.gra * obj.H);
            xg = obj.ChLength:-waveLen:0;
            drawGaugePoints@OpenChannel2d( obj, xg );
        end
    end
    
    methods( Access = protected, Static )
        function obtype = setOpenBoundaryCondition( )
            obtype = [ NdgEdgeType.Clamped, NdgEdgeType.SlipWall ];
        end
    end
    
end

