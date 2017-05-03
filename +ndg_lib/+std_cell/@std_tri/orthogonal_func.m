function [ fval ] = orthogonal_func(obj, N, ind, r, s, t)
%ORTHOGONALFUN 正交基函数
%   返回正交基函数在坐标（r,s）处函数值

% 坐标投影到矩阵
[a,b] = rstoab(r,s);

% 基函数序号转换为矩形内正交函数编号（i，j）
[i, j] = TransInd(N,ind);
fval = Simplex2DP(a,b,i,j);

end

