function [ mesh ] = read_mesh_file( N, casename, type )
%READ_MESH_FILE Summary of this function goes here
%   Detailed explanation goes here

% 读取二维几何网格
switch type
    case ndg_lib.std_cell_type.Tri
        cell = ndg_lib.std_cell.tri(N);
        mesh = ndg_lib.mesh.tri_mesh(cell, casename);
    case ndg_lib.std_cell_type.Quad
        cell = ndg_lib.std_cell.quad(N);
        mesh = ndg_lib.mesh.quad_mesh(cell, casename);
end% func

end


