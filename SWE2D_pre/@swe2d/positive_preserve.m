function [ f_Q ] = positive_preserve( obj, f_Q )
%POSITIVE_PRESERVE 修正水深与流量保证水深为正值
%   Detailed explanation goes here

hc = obj.mesh.cell_mean(f_Q(:,:,1));
qxc = obj.mesh.cell_mean(f_Q(:,:,2));
qyc = obj.mesh.cell_mean(f_Q(:,:,3));

[f_Q(:,:,1), f_Q(:,:,2), f_Q(:,:,3)] = ppreserve(obj.hmin, ...
    f_Q(:,:,1), f_Q(:,:,2), f_Q(:,:,3), hc, qxc, qyc);
end