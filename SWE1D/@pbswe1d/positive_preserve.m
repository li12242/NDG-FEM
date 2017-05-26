function [ f_Q ] = positive_preserve( obj, f_Q )
%POSE_PRESERVE Summary of this function goes here
%   Detailed explanation goes here

h = f_Q(:,:,1) - obj.bot;
hc = obj.mesh.cell_mean( h );
qc = obj.mesh.cell_mean( f_Q(:,:,2) );

[h, f_Q(:,:,2)] = ppreserve(obj.hmin, h, f_Q(:,:,2), hc, qc);
f_Q(:,:,1) = h + obj.bot;
end