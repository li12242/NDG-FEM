classdef mesh2d < gq_lib.mesh.mesh & ndg_lib.mesh.mesh2d
    %MESH2D Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = mesh2d(cell, Nv, vx, vy, K, EToV, EToR, EToBS)
            obj = obj@ndg_lib.mesh.mesh2d(cell, ...
                Nv, vx, vy, K, EToV, EToR, EToBS);
            vz = zeros(size(vx));
            
            obj = obj@gq_lib.mesh.mesh(cell, ...
                Nv, vx, vy, vz, K, EToV, EToR, EToBS);
        end
    end
    
    methods(Static)
        [Nv, vx, vy, K, EToV, EToR, EToBS] = read_from_file(casename);
    end
    
    methods(Hidden, Access=protected)
        function [nxq, nyq, nzq, wJsq] = face_scal(obj, vx, vy, vz, EToV)
            nxq = zeros(obj.cell.Nfqtotal, obj.K);
            nyq = zeros(obj.cell.Nfqtotal, obj.K);
            nzq = zeros(obj.cell.Nfqtotal, obj.K);
            % start index of each face node
            faceIndexStart = ones(obj.cell.Nface, 1); 
            for f = 2:obj.cell.Nface
                faceIndexStart(f) = faceIndexStart(f-1) + obj.cell.Nfq(f);
            end

            for f = 1:obj.cell.Nface
                Nfp = obj.cell.Nfq(f);
                face_x1 = vx(EToV(obj.cell.FToV(1,f), :))';
                face_x2 = vx(EToV(obj.cell.FToV(2,f), :))';
                face_y1 = vy(EToV(obj.cell.FToV(1,f), :))';
                face_y2 = vy(EToV(obj.cell.FToV(2,f), :))';

                ind = faceIndexStart(f):(faceIndexStart(f)+Nfp-1);
                nxq(ind, :) = repmat( (face_y2 - face_y1), Nfp, 1 );
                nyq(ind, :) = repmat(-(face_x2 - face_x1), Nfp, 1 );
            end

            % normalise
            Js = sqrt(nxq.*nxq+nyq.*nyq); 
            nxq = nxq./Js; 
            nyq = nyq./Js;
            wJsq = bsxfun(@times, obj.cell.wbq, Js.*0.5);
        end% func
    end
    
end

