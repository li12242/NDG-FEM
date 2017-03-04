function [ hlim ] = VB2d_JK( mesh, h )
%VB2D_JK Summary of this function goes here
%   Detailed explanation goes here

shape = mesh.Shape;
hlim  = Utilities.Limiter.Limiter2D.VB2d_JK_Mex...
            (h, mesh.J, shape.M, mesh.EToV, ...
            mesh.Shape.Fmask, mesh.Nv, mesh.x, mesh.y);
end

