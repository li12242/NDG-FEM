function [ hlim ] = VB2d_VA( mesh, h )
%VB2D_VA Summary of this function goes here
%   Detailed explanation goes here

shape = mesh.Shape;
hlim  = Utilities.Limiter.Limiter2D.VB2d_VA_Mex...
            (h, mesh.J, shape.M, mesh.EToV, ...
            mesh.Shape.Fmask, mesh.Nv, mesh.x, mesh.y);
end

