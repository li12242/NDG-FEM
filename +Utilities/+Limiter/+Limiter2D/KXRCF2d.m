function [ ind ] = KXRCF2d( mesh, h, u, v, distol )
%KXRCF2D Summary of this function goes here
%   Detailed explanation goes here

% distol = 1.0;
shape = mesh.Shape;
ind  = Utilities.Limiter.Limiter2D.KXRCF_detector2d_Mex...
    (h, u, v, mesh.J, mesh.sJ, ...
    shape.M, shape.Fmask, mesh.vmapM, mesh.vmapP, ...
    shape.Mes, mesh.x, mesh.y, mesh.nx, mesh.ny, distol);
end

