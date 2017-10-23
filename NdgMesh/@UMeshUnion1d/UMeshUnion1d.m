classdef UMeshUnion1d < UMeshUnion
    %LINE_MESH Summary of this class goes here
    %   Detailed explanation goes here

    methods(Hidden, Access = protected)
        function [rx, ry, rz, sx, sy, sz, tx, ty, tz, J] = assembleJacobiFactor(obj)
            xr = obj.cell.Dr*obj.x;
            J = xr; rx = 1./J;

            ry = zeros(size(rx));
            rz = zeros(size(rx)); 
            sx = zeros(size(rx));
            sy = zeros(size(rx));
            sz = zeros(size(rx)); 
            tx = zeros(size(rx));
            ty = zeros(size(ry));
            tz = zeros(size(rx)); 
        end
        function [nx, ny, nz, Js] = assembleFacialJaobiFactor( obj )
            xb = obj.x(obj.cell.Fmask, :);
            nx = ones(obj.cell.TNfp, obj.K);
            % Define outward normals
            [~, ind] = min(xb);
            nx(ind, :) = -1;
            Js = ones(size(nx));
            ny = zeros(obj.cell.Nface, obj.K);
            nz = zeros(obj.cell.Nface, obj.K);
        end
        
        function faceId = assembleGlobalFaceIndex( obj )
            faceId = zeros(obj.cell.Nface, obj.K);
            for f = 1:obj.cell.Nface
                faceId(f, :) = obj.EToV(obj.cell.FToV(1,f), :);
            end
        end% func
    end% methods
    
    methods
        obj = refine(obj, refine_level);
        
        function obj = UMeshUnion1d(cell, Nv, vx, K, EToV, EToR, BCToV)
            vy = zeros(size(vx)); % vy is all zeros
            vz = zeros(size(vx)); % vz is all zeros            
            obj = obj@UMeshUnion(cell, ...
                Nv, vx, vy, vz, K, EToV, EToR, BCToV);
        end% func
        
        function draw(obj)
            plot(obj.x, zeros(obj.cell.Np, obj.K), '.-');
        end
    end
    
end

