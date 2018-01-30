classdef NdgVertLimiter < NdgAbstractLimiter
    
    properties( SetAccess = protected )
        %> Number of vertices
        Nv
        %> number of cells connecting at each vertex
        Nvc
        %> maximum number of cells connecting at one vertex
        Nvcmax
        %> cell index containing each vertex
        VToK
        %> mesh index of cells containing each vertex
        VToM
        %> reconstruction weights of each cell connecting to the vertex
        VToW
    end
    
    methods
        function obj = NdgVertLimiter( mesh )
            obj = obj@NdgAbstractLimiter( mesh );
            obj.Nv = mesh(1).Nv;
            [ obj.Nvcmax, obj.Nvc, obj.VToK, obj.VToM, obj.VToW ] = ...
                assembleVertexCellConnect( obj );
        end
    end
    
    methods( Abstract )
        fphys = matLimit( obj, fphys, fieldId );
    end
    
    methods( Hidden)%, Access = protected )
        %> Calculate the averaged, maximum and minimum vertex value
        [ fvert, fvmin, fvmax, cvar ] = matEvaluateVertAverage( obj, fphys, fieldId );
        %> Find the cell
        [ Nvcmax, Nvc, VToK, VToM, VToW ] = assembleVertexCellConnect( obj )
    end
    
end

