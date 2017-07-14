function map = GetBoundMapList(obj,BC,phyID)
% Input:    
%   phyID:  physical ID
%   BC:     boundary condition matrix, size [nBC x 3]
%           elementaries: [physID, node1, node2]
%   

bclist = (BC(:,1)==phyID);
partialBC = BC(bclist, :);
nBC = size(partialBC,1); Nv = obj.Nv;
VToBF = spalloc(Nv, nBC, 2*nBC);
for ib = 1:nBC
    VToBF(partialBC(ib,[2,3]),ib) = 1;
end
FToBF = obj.SpFToV*VToBF;
[faces1, ~] = find(FToBF == 2);

Nfaces = obj.Shape.nFace;
element1 = floor( (faces1-1)/Nfaces )  + 1; 
face1    =   mod( (faces1-1), Nfaces ) + 1;

Nfp = obj.Shape.nFaceNode/obj.Shape.nFace; 
map = zeros(Nfp, nBC);
for i=1:nBC
    localfacelist = obj.Shape.GetFaceListAtFace(face1(i));
    facelist = (element1(i)-1)*obj.Shape.nFaceNode + localfacelist;
    map(:,i) = facelist;
end
end% function