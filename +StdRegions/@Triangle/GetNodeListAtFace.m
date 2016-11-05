function [ Nfp, nodelist ] = GetNodeListAtFace( obj, ind )
%GETNODELISTATFACE 获取第ind条边上节点全局节点编号与节点个数
%   节点按照逆时针排列

order = obj.nOrder;
Nfp = order+1;
facelist = obj.GetFaceListAtFace(ind);
faceListToNodeList = obj.GetFaceListToNodeList();
nodelist = faceListToNodeList(facelist);

end