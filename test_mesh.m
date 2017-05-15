%% line
line = ndg_lib.std_cell.line(2);

K = 100;
Nv = K+1;
EToV = [1:K; 2:(K+1)];
vx = linspace(0, 1, Nv)';
EToR = ones(K, 1);
EToBS = zeros(size(EToV)); 
EToBS([1, end]) = ndg_lib.bc_type.ZeroGrad;

mesh = ndg_lib.mesh.line_mesh(line, Nv, vx, K, EToV, EToR, EToBS);

%% triangle mesh
tri = ndg_lib.std_cell.tri(2);

EToV = [1,2,3; 3,2,4]';
vx = [0, 1, 0, 1]';
vy = [0, 0, 1, 1]';
K = 2;
Nv = 4;
EToR = [1, 1]';
EToBS = [5, 1, 5; 1, 5, 5]';

mesh = ndg_lib.mesh.tri_mesh(tri, Nv, vx, vy, K, EToV, EToR, EToBS);

casename = 'Conv2D/mesh/triangle';
mesh = ndg_lib.mesh.tri_mesh(tri, casename);

%% quad
quad = ndg_lib.std_cell.quad(4);

EToV = [1,2,5,4; 2,3,6,5; 4,5,8,7]';
vx = [0, 1, 2, 0, 1, 2, 0, 1]';
vy = [0, 0, 0, 1, 1, 1, 2, 2]';
K = 3;
Nv = 8;
EToR = [1, 1, 1]';
EToBS = [5,1,1,5; 1,5,1,1; 1,5,5,5]';

mesh = ndg_lib.mesh.quad_mesh(quad, Nv, vx, vy, K, EToV, EToR, EToBS)
