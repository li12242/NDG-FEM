%> @brief Here we have a brief description of the class.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDGOM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef NdgInnerEdge < handle
    
    properties( SetAccess = protected )
        %> mesh obj
        mesh
        %> edge std cell obj
        bcell
        %> num of edges
        Ne
        %> num of face nodes
        Nfp
        %> vertex index on each edge
        FToV
        %> local and adjacent cell index
        FToE
        %> local face index of local and adjacent cell
        FToF
        %> face to mesh index
        FToM
        %> edge type
        ftype
        
        %> interp node index of 1st ele on each edge
        FToN1
        %> interp node index of 2nd ele on each edge
        FToN2
        %> outward normal vector
        nx
        %> outward normal vector
        ny
        %> outward normal vector
        nz
        %> determination of edge Jacabian
        Js 
    end
    
    methods( Access = public )
        function obj = NdgInnerEdge( mesh, meshId )
            % set mesh object
            obj.mesh = mesh;
            
            % set mesh id
            obj.FToM = meshId;
            
            % set reference cell
            [ obj.bcell ] = obj.setEdgeReferCell( mesh );
            
            % set face node num
            obj.Nfp = obj.bcell.Np;
            
            % connect edge to elements
            [ obj.Ne, obj.FToE, obj.FToF, obj.FToV, obj.ftype ] ...
                = obj.assembleEdgeConnect( mesh );
            
            % connect node
            [ obj.FToN1, obj.FToN2, obj.nx, obj.ny, obj.nz, obj.Js ] ...
                = assembleNodeProject( obj, mesh );
        end
        %> evaluate R.H.S. for surface integral term
        frhs = matEvaluateStrongFromEdgeRHS( obj, fluxM, fluxP, fluxS );
    end
    
    methods( Abstract, Static, Access = protected )
        [ bcell ] = setEdgeReferCell( mesh );
    end
    
    methods( Access = private, Sealed )
        [ Nedge, FToE, FToF, FToV, ftype ] = assembleEdgeConnect( obj, mesh )
        [ FToN1, FToN2, nx, ny, nz, Js ] = assembleNodeProject( obj, mesh )
    end
    
end

