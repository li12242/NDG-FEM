function [rx, ry, rz, sx, sy, sz, tx, ty, tz, J] = ele_vol_factor( obj )
%ELE_VOL_FACTOR Summary of this function goes here
%   Detailed explanation goes here

xr = obj.cell.Dr*obj.x;
J = xr; rx = 1./J;

ry = zeros(size(rx));
rz = zeros(size(rx)); 
sx = zeros(size(rx));
sy = zeros(size(rx));
sz = zeros(size(rx)); 
tx = zeros(size(rx));
ty = zeros(size(ry));
tz = zeros(size(rx)); 
end

