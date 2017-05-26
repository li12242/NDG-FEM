function [ flux ] = flux_term( obj, f_Q )
%FLUX_TERM Summary of this function goes here
%   Detailed explanation goes here

flux = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);
[flux(:,:,1), flux(:,:,2)] = nodal_flux(obj.hmin, obj.gra, ...
    f_Q(:,:,1), f_Q(:,:,2), obj.bot, obj.mesh.EToR);
end

