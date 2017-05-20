function [ mesh ] = uniform_mesh( N, M, type )
%UNIFORM_MESH Summary of this function goes here
%   Detailed explanation goes here

xmin = 0; xmax = 1; 
ymin = 0; ymax = 1;
face_type = [ndg_lib.bc_type.Clamped,...
    ndg_lib.bc_type.Clamped, ...
    ndg_lib.bc_type.Clamped, ...
    ndg_lib.bc_type.Clamped];

switch type
    case ndg_lib.std_cell_type.Tri
        cell = ndg_lib.ndg_cell(N, type);
        [K,EToV,Nv,VX,VY,EToBS,EToR] = ...
            ndg_utility.uniform_mesh.tri_mesh(M, M, ...
            xmin, xmax, ymin, ymax, face_type);
        mesh = ndg_lib.mesh.tri_mesh(cell, Nv, VX, VY, K, EToV, EToR, EToBS);
        
    case ndg_lib.std_cell_type.Quad
        cell = ndg_lib.ndg_cell(N, type);
        [K,EToV,Nv,VX,VY,EToBS,EToR] = ...
            ndg_utility.uniform_mesh.quad_mesh(M, M, ...
            xmin, xmax, ymin, ymax, face_type);
        mesh = ndg_lib.mesh.quad_mesh(cell, Nv, VX, VY, K, EToV, EToR, EToBS);
end

end

