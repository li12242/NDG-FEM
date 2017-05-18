function [ dflux ] = hll_surf_term( obj, f_Q )
%HLL_SURF_TERM Summary of this function goes here
%   Detailed explanation goes here

dflux = zeros(obj.mesh.cell.Nfptotal, obj.mesh.K, obj.Nfield);
Nptol = ( obj.mesh.cell.Np*obj.mesh.K );

etaM = f_Q( obj.mesh.eidM );
etaP = f_Q( obj.mesh.eidP );
bot_M = obj.bot( obj.mesh.eidM ); q_M = f_Q( obj.mesh.eidM+Nptol );
bot_P = obj.bot( obj.mesh.eidP ); q_P = f_Q( obj.mesh.eidP+Nptol );

[dflux(:,:,1), dflux(:,:,2)] = hll_flux(obj.hmin, obj.gra, ...
    etaM, q_M, bot_M, etaP, q_P, bot_P, obj.mesh.nx);
end
