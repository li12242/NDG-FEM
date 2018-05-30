function [rx, sx, tx, ry, sy, ty, rz, sz, tz, J] = Jacobian3d(obj, cell)

xr = cell.Dr * obj.x; xs = cell.Ds * obj.x; xt = cell.Dt * obj.x;
yr = cell.Dr * obj.y; ys = cell.Ds * obj.y; yt = cell.Dt * obj.y;
zr = cell.Dr * obj.z; zs = cell.Ds * obj.z; zt = cell.Dt * obj.z;
J = xr .* (ys .* zt - zs .* yt ) - yr .* ( xs .* zt - zs .* xt ) ...
    + zr .* (xs .* yt - ys .* xt);

rx = ( ys .* zt - zs .* yt ) ./ J;
sx = - (yr .* zt - zr .* yt) ./ J;
tx = ( yr .* zs - zr .* ys ) ./ J;
ry = - ( xs .* zt - zs .* xt ) ./ J;
sy = ( xr .* zt - zr .* xt ) ./ J;
ty = - ( xr .* zs - zr .* xs ) ./ J;
rz = ( xs .* yt - ys .* xt ) ./ J;
sz = - ( xr .* yt - yr .* xt ) ./ J;
tz = ( xr .* ys - yr .* xs ) ./ J;
end