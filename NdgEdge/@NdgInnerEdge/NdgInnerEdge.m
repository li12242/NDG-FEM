%> @brief Here we have a brief description of the class.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDGOM software.
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef NdgInnerEdge < handle
    
    properties ( SetAccess = protected )
        %> mesh obj
        mesh
        %> mass matrix
        M
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
        nx, ny, nz
        %> determination of edge Jacabian
        Js
    end
    
    methods ( Access = public )
        function obj = NdgInnerEdge( meshUnion, meshId )
            obj.mesh = meshUnion(meshId); % set mesh object
            obj.FToM = meshId; % set mesh id
            [ obj.Nfp, obj.M ] = assembleMassMatrix( obj );
            % obj.eCell = obj.setEdgeReferCell( obj.mesh ); % set reference cell
            
            % connect edge to elements
            obj = obj.assembleEdgeConnect( obj.mesh );
            
            % connect node
            obj = assembleNodeProject( obj, obj.mesh );
        end
        
        %> evaluate R.H.S. for surface integral term
        frhs = matEvaluateStrongFromEdgeRHS( obj, fluxM, fluxP, fluxS );
        
        [ fM, fP ] = matEvaluateSurfValue( obj, fphys );
    end
    
    methods ( Abstract, Access = protected )
        %> connect edge to elements
        [ Nfp, M ] = assembleMassMatrix( obj );
        obj = assembleEdgeConnect( obj, mesh )
        obj = assembleNodeProject( obj, mesh )
    end
    
end

