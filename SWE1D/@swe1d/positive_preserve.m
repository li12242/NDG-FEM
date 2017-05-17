function [ f_Q ] = positive_preserve( obj, f_Q )
%POSE_PRESERVE Summary of this function goes here
%   Detailed explanation goes here

hc = obj.mesh.cell_mean(f_Q(:,:,1));
qc = obj.mesh.cell_mean(f_Q(:,:,2));

[f_Q(:,:,1), f_Q(:,:,2)] = ppreserve(obj.hmin, f_Q(:,:,1), f_Q(:,:,2), hc, qc);
end