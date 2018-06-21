function obj = Jacobian3d(obj, cell)

xr = cell.Dr * obj.x; xs = cell.Ds * obj.x; xt = cell.Dt * obj.x;
yr = cell.Dr * obj.y; ys = cell.Ds * obj.y; yt = cell.Dt * obj.y;
zr = cell.Dr * obj.z; zs = cell.Ds * obj.z; zt = cell.Dt * obj.z;
J = xr .* (ys .* zt - zs .* yt ) - yr .* ( xs .* zt - zs .* xt ) ...
    + zr .* (xs .* yt - ys .* xt);

obj.rx = ( ys .* zt - zs .* yt ) ./ J;
obj.sx = - (yr .* zt - zr .* yt) ./ J;
obj.tx = ( yr .* zs - zr .* ys ) ./ J;
obj.ry = - ( xs .* zt - zs .* xt ) ./ J;
obj.sy = ( xr .* zt - zr .* xt ) ./ J;
obj.ty = - ( xr .* zs - zr .* xs ) ./ J;
obj.rz = ( xs .* yt - ys .* xt ) ./ J;
obj.sz = - ( xr .* yt - yr .* xt ) ./ J;
obj.tz = ( xr .* ys - yr .* xs ) ./ J;

obj.J = J;
obj.Jz = zt;
end