classdef OpenChannel2d_TransposeBC < OpenChannel2d
    methods
        function obj = OpenChannel2d_TransposeBC( N, M, type )
            obj = obj@OpenChannel2d( N, M, type );
        end
        
        function drawGaugePoints( obj )
            xg = linspace( 0, obj.ChLength, 5 );
            drawGaugePoints@OpenChannel2d( obj, xg );
        end
    end
    
    methods( Access = protected, Static )
        function obtype = setOpenBoundaryCondition( )
            obtype = [ NdgEdgeType.ClampedDepth, NdgEdgeType.ClampedDepth ];
        end
    end
end