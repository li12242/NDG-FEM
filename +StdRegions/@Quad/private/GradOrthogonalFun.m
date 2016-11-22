function [ fdr, fds ] = GradOrthogonalFun( r,s,order,ind )
%GRADORTHOGONALFUN Summary of this function goes here
%   Detailed explanation goes here

% 转换编号为（i,j）形式
[i,j] = TransInd(order,ind);
% 计算正交基函数函数值

fdr = Polylib.GradJacobiP(r(:), 0, 0, i).*Polylib.JacobiP(s(:), 0, 0, j);
fds = Polylib.JacobiP(r(:), 0, 0, i).*Polylib.GradJacobiP(s(:), 0, 0, j);

end

