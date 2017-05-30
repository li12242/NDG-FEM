function [ mesh, obcfile ] = read_mesh_file( cell, casename )
%READ_MESH_FILE Summary of this function goes here
%   Detailed explanation goes here

% 读取二维几何网格
switch cell.type
    case ndg_lib.std_cell_type.Tri
        mesh = ndg_lib.mesh.tri_mesh(cell, casename);
    case ndg_lib.std_cell_type.Quad
        mesh = ndg_lib.mesh.quad_mesh(cell, casename);
end% func

% 生成开边界条件
hin = 1.0;
uin = 8.57;
qxin = hin*uin;
vert = read_edge_file(casename);
Nv = numel(vert);
time = 0;
Nt = numel(time);
Nfield = 3;
f_Q = zeros(Nv, Nfield, Nt);
f_Q(:, 1) = hin;
f_Q(:, 2) = qxin;
obcfile = ndg_lib.phys.obc_file();
obcfile.make_obc_file([casename, '.nc'], time, vert, f_Q);
obcfile.set_file([casename, '.nc']);
end

function vert = read_edge_file(casename)
filename = [casename, '.edge'];
fp = fopen(filename);
% read total number of 
Nf = fscanf(fp, '%d', 1);
fgetl(fp); % pass the rest of first line

% read face to vertex list
data = fscanf(fp, '%d %d %d %d', [4, Nf]);
ind = data([2,3], :);
vert = unique( ind(:) );

fclose(fp);
end% func

