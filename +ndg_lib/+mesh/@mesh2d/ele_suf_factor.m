function [nx, ny, nz, Js] = ele_suf_factor(obj, vx, vy, EToV)
%MESH_EDGE_FACTOR 计算单元内各个边外法线向量与Jacobian变换系数
%   计算单元边界应均为直线，每个边上各个节点外法线向量只有唯一值
%

nx = zeros(obj.cell.Nface, obj.K);
ny = zeros(obj.cell.Nface, obj.K);
nz = zeros(obj.cell.Nface, obj.K);

for f = 1:obj.cell.Nface
    face_x1 = vx(EToV(obj.cell.FToV(1,f), :));
    face_x2 = vx(EToV(obj.cell.FToV(2,f), :));
    face_y1 = vy(EToV(obj.cell.FToV(1,f), :));
    face_y2 = vy(EToV(obj.cell.FToV(2,f), :));
    
    nx(f, :) =  (face_y2 - face_y1);
    ny(f, :) = -(face_x2 - face_x1);
end

% normalise
Js = sqrt(nx.*nx+ny.*ny); 
nx = nx./Js; 
ny = ny./Js;
Js = Js.*0.5;
end

