classdef NdgSubMesh1d < NdgMesh1d
    
    properties( SetAccess = protected )
        nxInner
        nyInner
        nzInner
        flen
    end
    
    methods
        function obj = NdgSubMesh1d( cell, Nv, vx, K, EToV, EToR, BCToV )
            obj = obj@NdgMesh1d( cell, Nv, vx, K, EToV, EToR, BCToV );
            [ obj.nxInner, obj.nyInner, obj.nzInner, obj.flen ] = ...
                cell.matEvaluateInnerControlEdgeInfo( obj, x, y, z );
        end
    end
    
end

