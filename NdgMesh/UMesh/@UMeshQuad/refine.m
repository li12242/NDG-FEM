function obj = refine( obj, multi_rate )
%REFINE Summary of this function goes here
%   Detailed explanation goes here

for m = 1:multi_rate
    
    [newNv, newVX, newVY] = new_vertex(obj);
    [newK, newEToV, newEToBS, newEToR] = new_elemet_info(obj);
    
    obj = ndg_lib.mesh.quad_mesh(obj.cell, 'variable',...
        {newNv, newVX, newVY, newK, newEToV, newEToR, newEToBS});
end
end

function [newK, newEToV, newEToBS, newEToR] = new_elemet_info(obj)
newK = 4*obj.K;
newEToV = zeros(obj.cell.Nv, newK);
newEToBS = zeros(obj.cell.Nv, newK);
newEToR = int8(ones(newK, 1));
edge = ndg_lib.mesh.edge(obj); % 获得边对象
CNToE = (1:obj.K) + obj.Nv + edge.Nedge; % 单元形心节点编号

inner = ndg_lib.bc_type.Inner;
for k = 1:obj.K
    ksid = (k-1)*4;
    v = obj.EToV(:, k); % 单元顶点编号
    bc = obj.EToBS(:, k);
    for f1 = 1:obj.cell.Nface
        f2 = mod(f1-1 + 3, obj.cell.Nface)+1; % 前一个面编号
        n1 = find( ( (edge.kM-k).^2+(edge.fM-f1).^2 == 0 ) | ...
            ( (edge.kP-k).^2+(edge.fP-f1).^2 == 0 ) );
        n2 = find( ( (edge.kM-k).^2+(edge.fM-f2).^2 == 0 ) | ...
            ( (edge.kP-k).^2+(edge.fP-f2).^2 == 0 ) );
        
        newEToV(:, ksid+f1 ) = ...
            [v(f1), n1+obj.Nv, CNToE(k), n2+obj.Nv]';
        newEToBS(:, ksid+f1 ) = [bc(f1), inner, inner, bc(f2)];
        newEToR(ksid+f1 ) = obj.EToR(k);
    end
end
end

function [newNv, newVX, newVY] = new_vertex(obj)
    xf = obj.face_mean(obj.x); % 计算每个单元边中点
    yf = obj.face_mean(obj.y);

    edge = ndg_lib.mesh.edge(obj); % 获得边对象
    xe = zeros(edge.Nedge, 1);
    ye = zeros(edge.Nedge, 1);
    for n = 1:edge.Nedge
        xe(n) = xf(edge.fM(n), edge.kM(n));
        ye(n) = yf(edge.fM(n), edge.kM(n));
    end

    xc = obj.cell_mean(obj.x); % 计算单元形心
    yc = obj.cell_mean(obj.y);
    newVX = [obj.vx; xe; xc'];
    newVY = [obj.vy; ye; yc'];
    
    newNv = obj.Nv + edge.Nedge + obj.K;
end
