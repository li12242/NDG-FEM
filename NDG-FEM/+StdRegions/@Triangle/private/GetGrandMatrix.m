function [ Dr, Ds ] = GetGrandMatrix( N,r,s,V )
%GETGRANDMATRIX 计算节点基函数关于坐标（r,s）的导数矩阵
%   其中导数矩阵满足 
%   Dr_{ij} = \frac{\partial l_j}{\partial r} |_{(r = r_i, s = s_i)}

Np = (N+1)*(N+2)/2;
Vdr = zeros(Np, Np);
Vds = zeros(Np, Np);
for ind = 1:Np
    [Vdr(:,ind), Vds(:,ind)] = GradOrthogonalFun( r,s,N,ind );
end% for

Dr = Vdr/V; Ds = Vds/V;
end
