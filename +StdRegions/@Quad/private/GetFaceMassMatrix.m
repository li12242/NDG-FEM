function [ Mes ] = GetFaceMassMatrix( obj )
%GETFACEMASSMATRIX Summary of this function goes here
%   Detailed explanation goes here

order = obj.nOrder;
line = StdRegions.Line(order);
Mes = zeros(obj.nNode, obj.nFaceNode);

for f = 1:obj.nFace
    % iFaceList: the No. of ibth boundary local face node list. 
    faceind = obj.GetFaceListAtFace(f);
    [~,nodeind] = obj.GetNodeListAtFace(f);
    Mes(nodeind,faceind) = line.M;
end

end

