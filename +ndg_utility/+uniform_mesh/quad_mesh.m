function [K,EToV,Nv,VX,VY,EToBS,EToR] = quad_mesh(Mx, My, ...
    xmin, xmax, ymin, ymax, bc_type)
%QUAD_MESH 生成均匀四边形网格。
%   四边形网格计算域大小为 [xmin, xmax] x [ymin, ymax]，其中 x，y 每个轴上
%   单元个数分别为 Mx 与 My。所有四边形均为矩形。
%   默认节点编号从左下角开始，首先沿 x 轴进行循环。设置边界面类型分别为底部、上部、
%   左侧和右侧四个边界，数组 EToBS 中（每列）包括边界面上两个节点编号与边界类型。
% 
%   输入参数
%   Mx, My     - x，y 坐标轴上单元个数；
%   xmin,xmax  - x 坐标轴范围；
%   ymin,ymax  - y 坐标轴范围；
%   bctype - 底部、上部、左侧和右侧四个边界条件；
%
% Author: li12242 Tianjin University

%% Parameters
Nx = Mx + 1; % number of elements along x coordinate
Ny = My + 1;
K = Mx * My;
Nv = Nx * Ny;
EToR = int8( ones(K, 1) )*ndg_lib.mesh_type.Normal;

%% Define vectex
% The vertex is sorted along x coordinate. (x coordinate counts first)
VX   = linspace(xmin, xmax, Nx) ;
VY   = linspace(ymin, ymax, Ny)';
VX   = repmat(VX, 1, Ny) ;
VY   = repmat(VY, 1, Nx)';
VX   = VX(:);
VY   = VY(:);

%% Define EToV
% The element is conuting along x coordinate
EToV = zeros(4, Mx*My);
ind = 1;
for i = 1:My
    for j = 1:Mx
        % vertex index
        v1 = Nx*(i-1) + j;
        v2 = Nx*(i-1) + j + 1;
        v3 = Nx*i + j;
        v4 = Nx*i + j + 1;
        % Counterclockwise
        EToV(:, ind)=[v1, v2, v4, v3]';
        ind = ind +1;
    end% for
end% for

[ EToBS ] = ndg_utility.uniform_mesh.uniform_bc( Mx, My, EToV, bc_type );

end% func