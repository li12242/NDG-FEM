function [ nodelist ] = GetVertexNodeList( obj )
%GETVERTEXNODELIST 获取顶点所在节点序号
%   Detailed explanation goes here

Nfp = obj.nOrder+1; % No. of node on sigle face
nodelist = [1, Nfp, obj.nNode, obj.nNode-Nfp+1]';
end