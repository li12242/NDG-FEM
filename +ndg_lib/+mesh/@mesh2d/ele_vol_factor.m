function [rx, ry, rz, sx, sy, sz, tx, ty, tz, J] = ele_vol_factor(obj)
%MESH_VOLUME_FACTOR Summary of this function goes here
%   Detailed explanation goes here

xr = obj.cell.Dr*obj.x; xs = obj.cell.Ds*obj.x;
yr = obj.cell.Dr*obj.y; ys = obj.cell.Ds*obj.y;
J = -xs.*yr + xr.*ys;

rx = ys./J; sx =-yr./J; 
ry =-xs./J; sy = xr./J; 

tx = zeros(size(rx));
ty = zeros(size(ry));
rz = zeros(size(rx)); 
sz = zeros(size(rx)); 
tz = zeros(size(rx)); 
end

