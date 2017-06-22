function [ mesh ] = uniform_mesh( N, K )
%UNIFORM_MESH Summary of this function goes here
%   Detailed explanation goes here
line = ndg_lib.std_cell.line( N );

Nv = K+1;
EToV = [1:K; 2:(K+1)];

x1 = 0; x2 = 2; % º∆À„”Ú
vx = linspace(x1, x2, Nv)';
EToR  = zeros(K, 1);
EToBS = int8(ones(size(EToV)))*ndg_lib.bc_type.Inner;
EToBS(1) = ndg_lib.bc_type.Clamped;
EToBS(end) = ndg_lib.bc_type.ZeroGrad;

mesh = ndg_lib.mesh.line_mesh(line, Nv, vx, K, EToV, EToR, EToBS);

end

