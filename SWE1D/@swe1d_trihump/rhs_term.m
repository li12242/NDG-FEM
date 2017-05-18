function [ rhs ] = rhs_term( obj, f_Q )
%RHS_TERM Summary of this function goes here
%   Detailed explanation goes here

[ dflux ] = hll_surf_term( obj, f_Q );
[ flux ] = flux_term( obj, f_Q );
[ sb ] = topo_sour_term( obj, f_Q );
[ sf ] = fric_sour_term( obj, f_Q );

rhs(:,:,1) = -obj.mesh.rx.*( obj.mesh.cell.Dr*flux(:, :, 1) ) + ...
    obj.mesh.cell.LIFT*( obj.mesh.eidfscal.*dflux(:,:,1) ) + ...
    sb(:, :, 1) + sf(:, :, 1);

rhs(:,:,2) = -obj.mesh.rx.*( obj.mesh.cell.Dr*flux(:, :, 2) ) + ...
    obj.mesh.cell.LIFT*( obj.mesh.eidfscal.*dflux(:,:,2) ) + ...
    sb(:,:,2) + sf(:, :, 2);
end

