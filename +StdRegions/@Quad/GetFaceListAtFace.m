function [ facelist ] = GetFaceListAtFace( obj,ind )
%GETFACELISTATFACE Summary of this function goes here
%   Detailed explanation goes here

Nfp = obj.nOrder + 1; % 每边节点个数
contour = Nfp*(ind-1); % 前n条边节点总数
facelist = (contour+1):(contour+Nfp); % 第ind条边节点编号

end

