function [x, y, rx, sx, ry, sy, J] = GetEleGeometric( obj, VX, VY )
%GETELEGEOMETRIC ͶӰ�����㵥Ԫ�ļ�����Ϣ
%   ���� [VX, VY] - ���㵥Ԫ���㣨��ʱ�룩
%   ����ֵ 
%       [x, y] - ���㵥Ԫ�ڵ�����
%       [rx, sx, ry, sy] - ����任ϵ��
%       J - Jacobi�任ϵ��

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

