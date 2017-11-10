function [ mesh ] = uniform_mesh( N, cell_type, xlim, ylim, Mx, My, bc_type )
%UNIFORM_MESH Summary of this function goes here
%   Detailed explanation goes here

uniform_info = {xlim, ylim, Mx, My, bc_type};
switch cell_type
    case ndg_lib.std_cell_type.Tri
        cell = ndg_lib.std_cell.tri(N);
        mesh = ndg_lib.mesh.tri_mesh(cell, 'uniform', uniform_info);
        
    case ndg_lib.std_cell_type.Quad
        cell = ndg_lib.std_cell.quad(N);
        mesh = ndg_lib.mesh.quad_mesh(cell, 'uniform', uniform_info);
end
end

