function [ facelist ] = GetFaceListAtFace( obj,ind )
%GETFACELISTATFACE 返回第ind条边上节点局部编号
%   局部编号顺序只考虑边界上节点，

Nfp = obj.nOrder + 1; % 每边节点个数
contour = Nfp*(ind-1); % 前n条边节点总数
facelist = (contour+1):(contour+Nfp); % 第ind条边节点编号

end

