function [ dflux ] = surf_term( obj, f_Q )
%SURF_TERM Summary of this function goes here
%   Detailed explanation goes here

dflux = upwind_flux(f_Q, obj.f_extQ, obj.u, obj.v, ...
    obj.mesh.nx, obj.mesh.ny, obj.mesh.eidM, obj.mesh.eidP, ...
    obj.mesh.eidtype);

% cM = f_Q( obj.mesh.eidM ); cP = f_Q( obj.mesh.eidP ); 
% uM = obj.u( obj.mesh.eidM ); uP = obj.u( obj.mesh.eidP ); 
% vM = obj.v( obj.mesh.eidM ); vP = obj.v( obj.mesh.eidP );
% c_extM = obj.f_extQ( obj.mesh.eidM );
% 
% cP = nei_node_val(cM, cP, c_extM, obj.mesh.eidtype);
% 
% fluxS = zeros([obj.mesh.cell.Nfptotal, obj.mesh.K]); % numerical flux
% nx = obj.mesh.nx; 
% ny = obj.mesh.ny;
% 
% ind = ( (nx.*uM+ny.*vM) >= 0); 
% EM = nx.*uM.*cM + ny.*vM.*cM; fluxS( ind) = EM( ind);
% EM = nx.*uP.*cP + ny.*vP.*cP; fluxS(~ind) = EM(~ind);
% 
% % 外法向通量偏差，fn - fn*
% dflux = nx.*uM.*cM + ny.*vM.*cM - fluxS;
% end
% 
% function cP = nei_node_val(cM, cP, c_extM, ftype)
% %ADJ_NODE_VAL Summary of this function goes here
% %   Detailed explanation goes here
% 
% %% inner edge - default
% % ind = (ftype == ndg_lib.bc_type.Inner);
% % cP(ind) = cP(ind);
% 
% %% slip wall and non-slip wall
% ind = (ftype == ndg_lib.bc_type.SlipWall);
% cP(ind) = cM(ind);
% 
% ind = (ftype == ndg_lib.bc_type.SlipWall);
% cP(ind) = cM(ind);
% 
% %% zero gradient
% ind = (ftype == ndg_lib.bc_type.ZeroGrad);
% cP(ind) = cM(ind);
% 
% %% clamped
% ind = (ftype == ndg_lib.bc_type.Clamped);
% cP(ind) = 2*c_extM(ind) - cM(ind);
% 
% end
