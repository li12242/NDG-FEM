function [ nx, ny, sJ ] = GetFaceGeometric( obj, x, y )
%GETFACEGEOMETRIC Summary of this function goes here
%   Detailed explanation goes here

% get Face Normal vector & surface jacobi factor
% Input:    x - node coordinate, size [nNode, nElement]
%           y - node coordinate, size [nNode, nElement]
% Output:   nx - outward vector
%           ny - outward vector
%           sJ - face jacobi factor
K = size(x, 2);
xr = obj.Dr*x; yr = obj.Dr*y; xs = obj.Ds*x; ys = obj.Ds*y;
% interpolate geometric factors to face nodes
Fmask = obj.GetFaceListToNodeList();
fxr = xr(Fmask, :); fxs = xs(Fmask, :); 
fyr = yr(Fmask, :); fys = ys(Fmask, :);
% build normals
Nfp = obj.nFaceNode/4;
Nf = 4;
nx = zeros(Nf*Nfp, K); ny = zeros(Nf*Nfp, K);
fid1 = (1:Nfp)'; fid2 = (Nfp+1:2*Nfp)'; fid3 = (2*Nfp+1:3*Nfp)';
fid4 = (3*Nfp+1:4*Nfp)';
% face 1
nx(fid1, :) =  fyr(fid1, :); ny(fid1, :) = -fxr(fid1, :);
% face 2
nx(fid2, :) =  fys(fid2, :); ny(fid2, :) = -fxs(fid2, :);
% face 3
nx(fid3, :) = -fyr(fid3, :); ny(fid3, :) =  fxr(fid3, :);
% face 4
nx(fid4, :) = -fys(fid4, :); ny(fid4, :) =  fxs(fid4, :);
% normalise
sJ = sqrt(nx.*nx+ny.*ny); nx = nx./sJ; ny = ny./sJ;
end

