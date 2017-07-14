function [ detector ] = TVB_detector2d( mesh, h, factor )
%TVB_DETECTOR2D Summary of this function goes here
%   Detailed explanation goes here

shape  = mesh.Shape;
% factor = 0.2;
detector = Utilities.Limiter.Limiter2D.TVB_detector2d_Mex...
    (h, mesh.J,shape.M, shape.Fmask, mesh.EToE, factor);

end

