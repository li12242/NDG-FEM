function [ EToBS ] = uniform_bc( Mx, My, EToV, face_type )
%UNIFORM_BC 设置均匀网格上边界节点及边界条件。
%   在均匀网格中，默认节点编号从左下角开始，首先沿 x 轴进行循环。返回边界面从
%   底部、上部、左侧和右侧四个方向进行循环，数组 EToBS 中（每列）包括边界面上
%   两个节点编号与边界类型。
% Input:
%   Mx, My      - x、y 轴方向上边界面个数；
%   face_type   - 底部、上部、左侧和右侧边界条件；
% Output:
%   BSToV       - 每个边界面上顶点编号值；
% Usages:
% 

% 统计每个开边界顶点序号
BSToV = ones((Mx+My)*2, 3);
Nx = Mx + 1; Ny = My + 1; % x 与 y 坐标轴上顶点个数
Nv = Nx * Ny;
% 底部边界
st = 1; 
se = st + Mx - 1;
BSToV(st:se,[1,2]) = [1:(Nx-1); 2:Nx]';
BSToV(st:se, 3) = face_type(1); % 边界类型
% 上部边界
st = se + 1; 
se = st + Mx - 1;
vs = Nx*(Ny-1)+1; 
ve = Nx*Ny;
BSToV(st:se,[1,2]) = [vs:(ve-1); (vs+1):ve]';
BSToV(st:se, 3) = face_type(2); % 边界类型
% 左侧边界
st = se + 1; 
se = st + My - 1;
vs = 1; 
ve = Nx*(Ny-1)+1; 
vstrid = Nx;
BSToV(st:se,[1,2]) = [vs:vstrid:(ve-vstrid); (vs+vstrid):vstrid:ve]';
BSToV(st:se, 3) = face_type(3); % 边界类型
% 右侧边界
st = se + 1; 
se = st + My - 1;
vs = Nx; 
ve = Nx*Ny; 
vstrid = Nx;
BSToV(st:se,[1,2]) = [vs:vstrid:(ve-vstrid); (vs+vstrid):vstrid:ve]';
BSToV(st:se, 3) = face_type(4); % 边界类型

% 给每个边界面赋予编号
ind = min( BSToV(:, [1,2]), [], 2)*Nv + max( BSToV(:, [1,2]), [], 2);
ftype = BSToV(:, 3);

% 转换为单元边界条件 EToBS
[Nface, K] = size(EToV);
EToBS = int8(ones(Nface, K)).*ndg_lib.bc_type.Inner; % initialize
for k = 1:K
    for f = 1:Nface
        v1 = EToV(f, k);
        v2 = EToV(mod(f, Nface)+1, k);
        t = min(v1, v2)*Nv + max(v1, v2);
        tnd = find( abs(ind - t)<1e-10 );
        if isempty(tnd)
            continue
        else
            EToBS(f, k) = ftype(tnd);
        end
    end
end
end

