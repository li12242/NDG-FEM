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
        %> local interp node coodinate
        r, s, t
        %> mass matrix
        M
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
            obj.mesh = meshUnion( meshId ); % set mesh
            obj.FToM = meshId; % set mesh id
            obj = assembleMassMatrix( obj );
            
            % connect edge to elements
            obj = obj.assembleEdgeConnect( meshUnion );
            % connect node
            obj = obj.assembleNodeProject( meshUnion );
        end
        
        %> evaluate R.H.S. for surface integral term
        [ frhs ] = matEvaluateStrongFromEdgeRHS( obj, fluxM, fluxP, fluxS );
        [ frhs ] = matEvaluateStrongFormEdgeCentralRHS( obj, fluxM, fluxP );
        [ frhs ] = matEvaluateStrongFormEdgeAlterRHS( obj, fluxM, fluxP );
        [ fM, fP ] = matEvaluateSurfValue( obj, fphys );        
    end
    
    methods ( Access = public, Abstract )
        [ fnode ] = proj_vert2node( obj, fvert );
    end
    
    methods ( Abstract, Access = protected )
        %> connect edge to elements
        obj = assembleMassMatrix( obj );
        obj = assembleEdgeConnect( obj, mesh )
        obj = assembleNodeProject( obj, mesh )
    end
    
end

