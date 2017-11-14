%> @brief The abstract limiter class.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDGOM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef NdgAbstractLimiter < handle
    
    properties( SetAccess = protected )
        %> Number of mesh objects
        Nmesh
        %> mesh object array
        meshUnion
    end
    
    methods
        function obj = NdgAbstractLimiter(mesh)
            [ obj.Nmesh ] = numel( mesh );
            [ obj.meshUnion ] = mesh;
        end
    end
    
    methods( Abstract )
        %> employ the limiter function to reconstruct the result
        fphys = matLimit( obj, fphys, fieldId );
    end
    
end

