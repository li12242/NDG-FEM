function [dr, ds, dt] = derivative_orthogonal_func( obj, N, ind, r, s, t)
%GRADORTHOGONALFUN 标准三角形内正交函数导数值

% 坐标投影到矩阵
[a,b] = rstoab(r,s);

% 基函数序号转换为矩阵内函数编号（i，j）
[i, j] = TransInd(N,ind);

% 计算对（r,s）坐标的导数
[dr, ds] = GradSimplex2DP(a,b,i,j);
dt = zeros(size(dr));
end