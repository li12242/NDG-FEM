function [ mesh ] = read_mesh_file( N, casename, type )
%READ_MESH_FILE Summary of this function goes here
%   Detailed explanation goes here

% 读取二维几何网格
switch type
    case ndg_lib.std_cell_type.Tri
        cell = ndg_lib.std_cell.tri(N);
        mesh = ndg_lib.mesh.tri_mesh(cell, casename);
    case ndg_lib.std_cell_type.Quad
        cell = ndg_test.cell_test.quad_fullquad(N);
        mesh = ndg_test.mesh_test.mesh2d_fullquad(cell, casename);
end% func

end


