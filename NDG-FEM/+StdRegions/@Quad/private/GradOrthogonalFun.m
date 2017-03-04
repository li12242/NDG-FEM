function [ fdr, fds ] = GradOrthogonalFun( r,s,order,ind )
%GRADORTHOGONALFUN Summary of this function goes here
%   Detailed explanation goes here

% ת�����Ϊ��i,j����ʽ
[i,j] = TransInd(order,ind);
% ������������������ֵ

fdr = Polylib.GradJacobiP(r(:), 0, 0, i).*Polylib.JacobiP(s(:), 0, 0, j);
fds = Polylib.JacobiP(r(:), 0, 0, i).*Polylib.GradJacobiP(s(:), 0, 0, j);

end

