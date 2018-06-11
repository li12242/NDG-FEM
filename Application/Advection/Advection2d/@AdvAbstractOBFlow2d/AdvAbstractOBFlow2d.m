classdef AdvAbstractOBFlow2d < AdvAbstractConstFlow2d
    
    methods ( Access = public )
        function obj = AdvAbstractOBFlow2d()
            obj = obj@AdvAbstractConstFlow2d();
        end
        
        function [ fM, fP ] = matImposeBoundaryCondition( obj, edge, nx, ny, fM, fP, fext )
            ind = ( edge.ftype == enumBoundaryCondition.Clamped );
            fP(:, ind) = fext(:, ind);
        end
    end
    
    
    methods ( Access = protected )
        function matUpdateExternalField( obj, time, fphys )
            for m = 1:obj.Nmesh
                edge = obj.meshUnion(m).BoundaryEdge;
                obj.fext{m} = obj.getExtFunc( edge.xb, edge.yb, time );
            end
        end
    end
end