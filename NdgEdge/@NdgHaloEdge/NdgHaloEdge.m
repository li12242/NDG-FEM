%> @brief Edge information on unstructured meshes.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDG-FEM software.
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef NdgHaloEdge < handle
    
    properties(SetAccess=protected)
        %> mesh obj
        mesh
        %> edge std cell obj
        eCell
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
        %> face type of each edge
        ftype
    end
    
    methods(Abstract, Hidden=true, Access=protected)
        [ eCell ] = setEdgeReferCell( obj, mesh );
    end
    
    methods
        %> \brief constructor for Halo edge
        %> \param[in] meshUnion input mesh vector
        %> \param[in] locMeshId local mesh indices
        function obj = NdgHaloEdge( meshUnion, locMeshId )
            obj.mesh = meshUnion( locMeshId );
            
            % set reference element
            obj.eCell = obj.setEdgeReferCell( obj.mesh );
            
            % set face node num
            obj.Nfp = obj.eCell.Np;
            
            % connect edge to elements
            [ obj.Ne, obj.FToE, obj.FToF, obj.FToV, obj.FToM, obj.ftype ] ...
                = obj.assembleEdgeConnect( obj.mesh, locMeshId );
            
            [ obj.FToN1, obj.FToN2, obj.nx, obj.ny, obj.nz, obj.Js ] ...
                = assembleNodeProject( obj, meshUnion );
        end
        
        
    end
    
    methods( Abstract, Access = protected )
        [ Nedge, FToE, FToF, FToV, FToM, ftype ] = assembleEdgeConnect( obj, mesh, meshId );
        [ FToN1, FToN2, nx, ny, nz, Js ] = assembleNodeProject( obj, meshUnion )
    end
end

