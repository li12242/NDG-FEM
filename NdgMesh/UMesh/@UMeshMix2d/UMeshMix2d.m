%> @brief Here we have a brief description of the class.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDG-FEM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef UMeshMix2d < handle
    
    properties
        triMesh
        quadMesh
    end
    
    methods
        function obj = UmeshMix2d(N, Nv, vx, vy, ...
                Ktri, TriToV, TriToR, Kquad, QuadToR, QuadToV, BCToV)
            tri = StdTri(N);
            quad = StdQuad(N);
            
            obj.triMesh = UMeshTri(tri);
        end% func
    end
    
end

