function [ rhsQ ] = rhs_term(obj, f_Q )
%RHS_TERM Summary of this function goes here
%   Detailed explanation goes here

% [ rhsQ ] = rhs_term(f_Q, obj.f_extQ, obj.u, obj.v, ...
%     obj.mesh.nx, obj.mesh.ny, ...
%     obj.mesh.eidM, obj.mesh.eidP, obj.mesh.eidtype, obj.mesh.EToR, ... 
%     obj.mesh.cell.Dr, obj.mesh.cell.Ds, ...
%     obj.mesh.rx, obj.mesh.ry, obj.mesh.sx, obj.mesh.sy, ...
%     obj.mesh.cell.LIFT, obj.mesh.J, obj.mesh.Js);

[ rhsQ ] = rhs_fv_term(f_Q, obj.f_extQ, obj.u, obj.v, ...
    obj.loc_fv.subcell.Nedge, ...
    obj.loc_fv.subcell.v1, obj.loc_fv.subcell.v2, ...
    obj.loc_fv.nx, obj.loc_fv.ny, obj.loc_fv.ds, ...
    obj.mesh.eidM, obj.mesh.eidP, ...
    obj.mesh.nx, obj.mesh.ny, ...
    obj.loc_fv.Js, obj.loc_fv.subcell.LIFT, ...
    obj.loc_fv.vol, obj.mesh.eidtype, obj.mesh.EToR);

% rhsQ = rhsQ + rhsV;


% [ rhsV ] = rhs_fv_term(f_Q, obj.u, obj.v, ...
%     obj.loc_fv.Nedge, obj.loc_fv.v1, obj.loc_fv.v2, ...
%     obj.loc_fv.nx, obj.loc_fv.ny, obj.loc_fv.ds, ...
%     obj.mesh.eidM, obj.mesh.eidP, ...
%     obj.mesh.nx, obj.mesh.ny, ...
%     obj.loc_fv.Js, obj.loc_fv.LIFT, obj.loc_fv.vol);
% 
% ind = (obj.mesh.EToR == ndg_lib.mesh_type.Refine);
% rhsQ(:, ind) = rhsV(:, ind);
% 

end


% function [ rhsV ] = rhs_fv_term(v_Q, u, v, ...
%     Nedge, v1, v2, nx, ny, ds, ...
%     eidM, eidP, bnx, bny, Js, LIFT, vol)
% % innter flux
% [ rhs1 ] = inner_flux(v_Q, u, v, Nedge, v1, v2, nx, ny, ds);
% [ rhs2 ] = bound_flux(v_Q, u, v, eidM, eidP, bnx, bny, Js, LIFT);
% 
% rhsV = (rhs1 + rhs2)./vol;
% 
% end% func
% 
% function [ rhs ] = inner_flux(v_Q, u, v, Nedge, v1, v2, nx, ny, ds)
% [Np, K] = size(v_Q);
% rhs = zeros(Np, K);
% for k = 1:K
%     for n = 1:Nedge
%         n1 = v1(n, k);
%         n2 = v2(n, k);
%         flux = upwind_flux( v_Q(n1), v_Q(n2), ...
%             u(n1), v(n1), nx(n, k), ny(n, k) );
%         rhs(n1) = rhs(n1) - flux * ds(n, k);
%         rhs(n2) = rhs(n2) + flux * ds(n, k);
%     end
% end
% end% func
% 
% function [ rhs ] = bound_flux(v_Q, u, v, eidM, eidP, nxM, nyM, Js, LIFT)
% 
% cM = v_Q(eidM); cP = v_Q(eidP);
% uM = u(eidM); vM = v(eidM);
% flux = upwind_flux(cM, cP, uM, vM, nxM, nyM);
% rhs = - LIFT*(Js.*flux);
% 
% end% func
% 
% function [ flux ] = upwind_flux(c1, c2, u, v, nx, ny)
% flux = zeros(size(c1));
% 
% ind = ( (u.*nx + v.*ny) >= 0 );
% flux(ind) = c1(ind).*( u(ind).*nx(ind) + v(ind).*ny(ind) );
% 
% ind = ( (u.*nx + v.*ny) < 0 );
% flux(ind) = c2(ind).*( u(ind).*nx(ind) + v(ind).*ny(ind) );
% end% func
