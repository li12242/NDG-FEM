function [ nodelist ] = GetVertexNodeList( obj )
%GETVERTEXNODELIST Summary of this function goes here
%   Detailed explanation goes here

Nfp = obj.nOrder+1;   % No. of node on sigle face
nodelist = [1, Nfp, obj.nNode]';
end

