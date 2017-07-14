function [obj] = Connect2D(obj,EToV)

Nfaces = obj.Shape.nFace;
% Find number of elements and vertices
K = obj.nElement; Nv = obj.Nv;
% Create face to node connectivity matrix
TotalFaces = Nfaces*K;
% List of local face to local vertex connections
vn = zeros(Nfaces, 2);
for f = 1:Nfaces
    vn(f, :) = [f, mod(f, Nfaces)+1];
end
% Build global face to node sparse array
obj.SpFToV = spalloc(TotalFaces, Nv, 2*TotalFaces);
sk = 1;
for k=1:K
  for face=1:Nfaces
    obj.SpFToV(sk, EToV(k, vn(face,:))) = 1;
    sk = sk+1;
  end
end
% Build global face to global face sparse array
SpFToF = obj.SpFToV*obj.SpFToV' - 2*speye(TotalFaces);
% Find complete face to face connections
[faces1, faces2] = find(SpFToF==2);
% Convert face global number to element and face numbers
element1 = floor( (faces1-1)/Nfaces ) + 1;
face1    =   mod( (faces1-1),Nfaces ) + 1;
element2 = floor( (faces2-1)/Nfaces ) + 1;
face2    =   mod( (faces2-1),Nfaces ) + 1;
% Rearrange into nElements x Nfaces sized arrays
ind = sub2ind([K, Nfaces], element1, face1);
EToE = (1:K)'*ones(1,Nfaces); EToF = ones(K,1)*(1:Nfaces);
EToE(ind) = element2; EToF(ind) = face2;
% Assignment
obj.EToE = EToE;
obj.EToF = EToF;

end% function