function [ facelist ] = GetFaceListAtFace( obj,ind )
%GETFACELISTATFACE Summary of this function goes here
%   Detailed explanation goes here

Nfp = obj.nOrder + 1; % ÿ�߽ڵ����
contour = Nfp*(ind-1); % ǰn���߽ڵ�����
facelist = (contour+1):(contour+Nfp); % ��ind���߽ڵ���

end

