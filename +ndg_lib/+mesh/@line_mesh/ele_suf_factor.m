function [nx, ny, nz, Js] = ele_suf_factor(obj, vx, vy, EToV)
%ELE_SUF_FACTOR Summary of this function goes here
%   Detailed explanation goes here

xb = obj.x(obj.cell.Fmask, :);
nx = ones(obj.cell.Nfptotal, obj.K);
% Define outward normals
[~, ind] = min(xb);
nx(ind, :) = -1;
Js = ones(size(nx));

ny = zeros(obj.cell.Nface, obj.K);
nz = zeros(obj.cell.Nface, obj.K);
end

