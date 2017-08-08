function [ fdr, fds, fdt ] = derivative_orthogonal_func(obj, N, ind, r, s, t)
%GRADORTHOGONALFUN Summary of this function goes here
%   Detailed explanation goes here

% 转换编号为（i,j）形式
[i,j] = trans_ind(N,ind);
% 计算正交基函数函数值

fdr = Polylib.GradJacobiP(r(:), 0, 0, i).*Polylib.JacobiP(s(:), 0, 0, j);
fds = Polylib.JacobiP(r(:), 0, 0, i).*Polylib.GradJacobiP(s(:), 0, 0, j);
fdt = zeros(size(fdr));
end

