function [ Nfp, nodelist ] = GetNodeListAtFace( obj, ind )
%GETNODELISTATFACE Summary of this function goes here
%   Detailed explanation goes here

Nfp = obj.nOrder+1;

facelist = obj.GetFaceListAtFace(ind);
list = obj.GetFaceListToNodeList();
nodelist = list(facelist);
end

