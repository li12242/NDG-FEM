classdef OpenChannel2d_TransposeBC < OpenChannel2d
    methods
        function obj = OpenChannel2d_TransposeBC( N, M )
            obj = obj@OpenChannel2d( N, M );
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
        
        %> Set exact solution for propagation wave
        function [h, hu] = setExactSolution( obj, x, time )
            [ h, hu ] = obj.setOBC( x, time );
        end% func
    end
end