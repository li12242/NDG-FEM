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
        %> face node index of local element
        FToFN1
        %> face node index of adjacent element
        FToFN2
        %> outward normal vector
        nx
        %> outward normal vector
        ny
        %> outward normal vector
        nz
        %> determination of edge Jacabian
        Js
        %> lift matrix
        LIFT
    end
    
    methods( Access = public )
        function obj = NdgInnerEdge( meshUnion, meshId )
            obj.mesh = meshUnion(meshId); % set mesh object
            obj.FToM = meshId; % set mesh id
            obj.eCell = obj.setEdgeReferCell( obj.mesh ); % set reference cell
            obj.Nfp = obj.eCell.Np; % set face node num
            obj.LIFT = obj.mesh.cell.LIFT;
            
            % connect edge to elements
            [ obj.Ne, obj.FToE, obj.FToF, obj.FToV ] ...
                = obj.assembleEdgeConnect( obj.mesh );
            
            % connect node
            [ obj.FToN1, obj.FToN2, obj.FToFN1, obj.FToFN2, ...
                obj.nx, obj.ny, obj.nz, obj.Js ] ...
                = assembleNodeProject( obj, obj.mesh );
        end
        
        %> evaluate R.H.S. for surface integral term
        frhs = matEvaluateStrongFromEdgeRHS( obj, fluxM, fluxP, fluxS );
        
        [ fM, fP ] = matEvaluateSurfValue( obj, fphys );
    end
    
    methods( Abstract, Static, Access = protected )
        %> set reference boundary element
        [ bcell ] = setEdgeReferCell( mesh );
    end
    methods( Abstract, Access = protected )
        %> connect edge to elements
        [ Nedge, FToE, FToF, FToV ] = assembleEdgeConnect( obj, mesh )
        [ FToN1, FToN2, FToFN1, FToFN2, nx, ny, nz, Js ] = assembleNodeProject( obj, mesh )
    end
    
end

