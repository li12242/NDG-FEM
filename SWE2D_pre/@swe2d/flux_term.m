function [ E, G ] = flux_term( obj, f_Q )
%FLUX_TERM Calculate the flux term for the volume calculation.
%   Detailed explanation goes here

E = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);
G = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);

[E(:,:,1), E(:,:,2), E(:,:,3), G(:,:,1), G(:,:,2), G(:,:,3)] = ...
    nodal_flux(obj.hmin, obj.gra, ...
    f_Q(:,:,1), f_Q(:,:,2), f_Q(:,:,3), ...
    obj.bot, obj.mesh.EToR);
end

