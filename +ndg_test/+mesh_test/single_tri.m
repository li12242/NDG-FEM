function [ mesh ] = single_tri( cell )
%SINGLE_TRI Summary of this function goes here
%   Detailed explanation goes here

if (cell.type ~= ndg_lib.std_cell_type.Tri)
    error('The input cell type is not triangle.');
end

EToV = [1,2,3]';
vx = [0, .55, 0.15]'/0.55*pi/2;
vy = [0, 0.15, .55]'/0.55*pi/2;
K = 1;
Nv = 3;
EToR = 1;
EToBS = [5, 1, 5]';

mesh = ndg_lib.mesh.tri_mesh(cell, ...
    Nv, vx, vy, ...
    K, EToV, EToR, EToBS);

end

