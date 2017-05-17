function [ mesh ] = uniform_mesh( N, K )
%MESH_GEN Summary of this function goes here
%   Detailed explanation goes here

line = ndg_lib.std_cell.line(N);

Nv = K+1; % 节点个数
EToV = [1:K; 2:(K+1)]; % 单元节点编号

x1 = 0; x2 = 1e3; % 计算域
vx = linspace(x1, x2, Nv)';
EToR  = zeros(K, 1);
EToBS = uint8(ones(size(EToV)))*ndg_lib.bc_type.Inner;
EToBS([1, end]) = ndg_lib.bc_type.ZeroGrad;

mesh = ndg_lib.mesh.line_mesh(line, Nv, vx, K, EToV, EToR, EToBS);
end

