function [ rhs ] = rhs_term( obj, f_Q )
%RHS_TERM 一维浅水方程强形式离散右端项计算
%   右端项包括体积积分 flux，面积分上法向通量之差 dflux，底坡源项 sb，
%   摩阻源项 sf 等。

[ dflux ] = hll_surf_term( obj, f_Q );
[ flux ] = flux_term( obj, f_Q );
[ sb ] = topo_sour_term( obj, f_Q );
[ sf ] = fric_sour_term( obj, f_Q );

rhs(:,:,1) = -obj.mesh.rx.*( obj.mesh.cell.Dr*flux(:, :, 1) ) + ...
    obj.mesh.cell.LIFT*( obj.mesh.Js.*dflux(:,:,1) )./obj.mesh.J + ...
    sb(:, :, 1) + sf(:, :, 1);

rhs(:,:,2) = -obj.mesh.rx.*( obj.mesh.cell.Dr*flux(:, :, 2) ) + ...
    obj.mesh.cell.LIFT*( obj.mesh.Js.*dflux(:,:,2) )./obj.mesh.J + ...
    sb(:,:,2) + sf(:, :, 2);
end

