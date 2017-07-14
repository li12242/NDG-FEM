function [nx, ny, nz, Js] = ele_suf_factor(obj, vx, vy, vz, EToV)
%MESH_EDGE_FACTOR 计算单元内各个边外法线向量与Jacobian变换系数
%   计算单元边界应均为直线，每个边上各个节点外法线向量只有唯一值
%

nx = zeros(obj.cell.Nfptotal, obj.K);
ny = zeros(obj.cell.Nfptotal, obj.K);
nz = zeros(obj.cell.Nfptotal, obj.K);

faceIndexStart = ones(obj.cell.Nface, 1); % start index of each face node
for f = 2:obj.cell.Nface
    faceIndexStart(f) = faceIndexStart(f-1) + obj.cell.Nfp(f-1);
end
            
for f = 1:obj.cell.Nface
    Nfp = obj.cell.Nfp(f);
    face_x1 = vx(EToV(obj.cell.FToV(1,f), :))';
    face_x2 = vx(EToV(obj.cell.FToV(2,f), :))';
    face_y1 = vy(EToV(obj.cell.FToV(1,f), :))';
    face_y2 = vy(EToV(obj.cell.FToV(2,f), :))';
    
    ind = faceIndexStart(f):(faceIndexStart(f)+Nfp-1);
    nx(ind, :) = repmat( (face_y2 - face_y1), Nfp, 1 );
    ny(ind, :) = repmat(-(face_x2 - face_x1), Nfp, 1 );
end

% normalise
Js = sqrt(nx.*nx+ny.*ny); 
nx = nx./Js; 
ny = ny./Js;
Js = Js.*0.5;
end

