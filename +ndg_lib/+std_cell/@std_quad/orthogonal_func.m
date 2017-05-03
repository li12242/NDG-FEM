function [ fval ] = orthogonal_func(obj, N, ind, r, s, t)
%ORTHOGONALFUN 计算四边形内正交基函数在节点处函数值
%   四边形内正交基函数通过一维正交函数展开而成

% 转换编号为（i,j）形式
[i,j] = TransInd(N,ind);
% 计算正交基函数函数值
fval = Polylib.JacobiP(r(:), 0, 0, i).*Polylib.JacobiP(s(:), 0, 0, j);

end

