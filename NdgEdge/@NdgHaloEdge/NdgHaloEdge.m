%> @brief Edge information on unstructured meshes.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDG-FEM software.
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef NdgHaloEdge < NdgInnerEdge
    
    properties ( SetAccess = protected )
        %> face type of each edge
        ftype
        %> boundary coordinate
        xb, yb, zb
    end
    
    methods
        %> \brief constructor for Halo edge
        %> \param[in] meshUnion input mesh vector
        %> \param[in] locMeshId local mesh indices
        function obj = NdgHaloEdge( meshUnion, locMeshId, BCToV )
            obj = obj@NdgInnerEdge( meshUnion, locMeshId );
            obj = assembleBoundaryConnection(obj, BCToV);
        end
        
        %> evaluate right-hand-side for surface integral term
        frhs = matEvaluateStrongFromEdgeRHS( obj, fluxM, fluxP, fluxS );
        %> get surface values from physical field
        [ fM, fP ] = matEvaluateSurfValue( obj, fphys );
        
    end
    
    methods( Abstract, Access = protected )
        obj = assembleBoundaryConnection(obj, BCToV);
    end
end

