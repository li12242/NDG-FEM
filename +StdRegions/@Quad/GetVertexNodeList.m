function [ nodelist ] = GetVertexNodeList( obj )
%GETVERTEXNODELIST ��ȡ�������ڽڵ����
%   Detailed explanation goes here

Nfp = obj.nOrder+1; % No. of node on sigle face
nodelist = [1, Nfp, obj.nNode, obj.nNode-Nfp+1]';
end