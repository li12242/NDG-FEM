function [ fdr, fds ] = GradOrthogonalFun( r,s,order,ind )
%GRADORTHOGONALFUN 标准三角形内正交函数导数值

% 坐标投影到矩阵
[a,b] = rstoab(r,s);

% 基函数序号转换为矩阵内函数编号（i，j）
[i, j] = TransInd(order,ind);

% 计算对（r,s）坐标的导数
[fdr, fds] = GradSimplex2DP(a,b,i,j);
end