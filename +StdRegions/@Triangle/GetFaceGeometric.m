function [nx, ny, sJ] = GetFaceGeometric(obj, x, y)
%GETFACEGEOMETRIC 计算投影到计算单元上边界几何信息
%   参数: [x,y] - 计算单元内节点坐标
%   返回值: [nx,ny] - 边界节点上外法线向量
%           sJ - 边界积分Jacobi变换系数（边界长度）

K = size(x, 2);
xr = obj.Dr*x; yr = obj.Dr*y; 
xs = obj.Ds*x; ys = obj.Ds*y;
% interpolate geometric factors to face nodes
Fmask = obj.GetFaceListToNodeList();
fxr = xr(Fmask, :); fxs = xs(Fmask, :); 
fyr = yr(Fmask, :); fys = ys(Fmask, :);
% build normals
Nfp = obj.nOrder+1;
nx = zeros(obj.nFaceNode, K); 
ny = zeros(obj.nFaceNode, K);
fid1 = (1:Nfp)'; 
fid2 = (Nfp+1:2*Nfp)'; 
fid3 = (2*Nfp+1:3*Nfp)';
% face 1
nx(fid1, :) =  fyr(fid1, :); 
ny(fid1, :) = -fxr(fid1, :);
% face 2
nx(fid2, :) =  fys(fid2, :)-fyr(fid2, :); 
ny(fid2, :) = -fxs(fid2, :)+fxr(fid2, :);
% face 3
nx(fid3, :) = -fys(fid3, :); 
ny(fid3, :) =  fxs(fid3, :);
% normalise
sJ = sqrt(nx.*nx+ny.*ny); 
nx = nx./sJ; 
ny = ny./sJ;

end

