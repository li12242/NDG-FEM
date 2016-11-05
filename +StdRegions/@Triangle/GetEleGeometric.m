function [x, y, rx, sx, ry, sy, J] = GetEleGeometric( obj, VX, VY )
%GETELEGEOMETRIC 投影到计算单元的几何信息
%   参数 [VX, VY] - 计算单元顶点（逆时针）
%   返回值 
%       [x, y] - 计算单元节点坐标
%       [rx, sx, ry, sy] - 坐标变换系数
%       J - Jacobi变换系数

assert((size(VX,1)==3), 'Tri GetEleGeometric: input vertex VX has fault')
assert((size(VY,1)==3), 'Tri GetEleGeometric: input vertex VY has fault')

x = 0.5*(-(obj.r+obj.s)*VX(1,:)+(1+obj.r)*VX(2,:)+(1+obj.s)*VX(3,:));
y = 0.5*(-(obj.r+obj.s)*VY(1,:)+(1+obj.r)*VY(2,:)+(1+obj.s)*VY(3,:));

xr = obj.Dr*x; xs = obj.Ds*x; 
yr = obj.Dr*y; ys = obj.Ds*y; 
J = -xs.*yr + xr.*ys;

rx = ys./J; sx =-yr./J; 
ry =-xs./J; sy = xr./J;

end

