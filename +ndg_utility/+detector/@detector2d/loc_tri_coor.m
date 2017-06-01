function [ rd, sd ] = loc_tri_coor( obj, kd, xd, yd )
%LOC_TRI_COOR Summary of this function goes here
%   Detailed explanation goes here

rd = ones(obj.Nd, 1);
sd = ones(obj.Nd, 1);
% 计算所在单元内三个面积坐标 L1、L2、L3
fmask1 = obj.mesh.cell.Fmask(1, :)';
fmask2 = obj.mesh.cell.Fmask(end, :)';
for n = 1:obj.Nd
    xv1 = obj.mesh.x( fmask1, kd(n) );
    yv1 = obj.mesh.y( fmask1, kd(n) );
    xv2 = obj.mesh.x( fmask2, kd(n) );
    yv2 = obj.mesh.y( fmask2, kd(n) );
    b = sqrt( (xv1 - xv2).^2 + (yv1 - yv2).^2 );
    a1 = sqrt( (xd(n) - xv1).^2 + (yd(n) - yv1).^2 );
    a2 = sqrt( (xd(n) - xv2).^2 + (yd(n) - yv2).^2 );
    d = ( b + a1 + a2 )./2;
    L = sqrt( d.*(d-b).*(d-a1).*(d-a2) )./obj.mesh.vol( kd(n) );
    rd(n) = -L(2) + L(3) - L(1);
    sd(n) = -L(2) - L(3) + L(1);
end

end
