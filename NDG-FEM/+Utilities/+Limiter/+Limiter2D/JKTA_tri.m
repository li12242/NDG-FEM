function hlim = JKTA_tri(mesh, h)
% Slope limiting from Jawahear and Kamath and modified by Tu and Aliabadi (2005)

% 1. compute geometric information for 4 element patch containing each element
% Build average matrix
AVE = sum(mesh.Shape.M)/2;

% Compute displacements from center of nodes for Taylor expansion of limited fields
dropAVE = eye(mesh.Shape.nNode)-ones(mesh.Shape.nNode,1)*AVE;
dx   = dropAVE*mesh.x; dy = dropAVE*mesh.y;

% Extract coordinates of vertices and centers of elements
ids = mesh.Shape.getVertexNodeList;
xv = mesh.x(ids, :);
yv = mesh.y(ids, :);

% Find neighbors in patch
E1 = mesh.EToE(:,1)'; E2 = mesh.EToE(:,2)'; E3 = mesh.EToE(:,3)';

% compute element centers
xc0 = AVE*mesh.x; xc1 = xc0(E1); xc2 = xc0(E2); xc3 = xc0(E3);
yc0 = AVE*mesh.y; yc1 = yc0(E1); yc2 = yc0(E2); yc3 = yc0(E3);

% compute face unit normals and lengths
fnx = [yv(2,:)-yv(1,:);yv(3,:)-yv(2,:);yv(1,:)-yv(3,:)]; 
fny = -[xv(2,:)-xv(1,:);xv(3,:)-xv(2,:);xv(1,:)-xv(3,:)];
fL = sqrt(fnx.^2 + fny.^2); fnx = fnx./fL; fny = fny./fL;

% Find boundary faces for each face 
BC = mesh.EToE - [1:mesh.nElement]'*ones(1, mesh.Shape.nFace);
id1 = find(~BC(:,1)); id2 = find(~BC(:,2)); id3 = find(~BC(:,3));


% Compute weights for face gradients 
A0 = AVE*mesh.J*2/3; A1 = A0+A0(E1); A2 = A0+A0(E2); A3 = A0+A0(E3);

% Compute location of centers of reflected ghost elements at boundary faces
H1 = 2*(A0(id1)./fL(1,id1)); 
xc1(id1) = xc1(id1) + 2*fnx(1,id1).*H1; yc1(id1) = yc1(id1) + 2*fny(1,id1).*H1; 

H2 = 2*(A0(id2)./fL(2,id2)); 
xc2(id2) = xc2(id2) + 2*fnx(2,id2).*H2; yc2(id2) = yc2(id2) + 2*fny(2,id2).*H2; 

H3 = 2*(A0(id3)./fL(3,id3)); 
xc3(id3) = xc3(id3) + 2*fnx(3,id3).*H3; yc3(id3) = yc3(id3) + 2*fny(3,id3).*H3; 

% Compute cell averages of conserved variables
VC0 = AVE*h;

% Find neighbor values of conserved variables
hC = VC0(mesh.EToE');

% find value of primitive variables in patches
VC1 = hC(1,:);  VC2 = hC(2,:); VC3 = hC(3,:);
Nfp = mesh.Shape.nFaceNode/mesh.Shape.nFace;
ids = [1;Nfp;Nfp+1;2*Nfp;2*Nfp+1;3*Nfp]; % facelist of 3 vertex
VA = (h(mesh.vmapM(ids, :)) + h(mesh.vmapP(ids, :)) )/2;

% Compute values of centers of reflected ghost elements at boundary faces
VC1(id1) = (VA(1,id1) + VA(2,id1)) - VC0(id1);
VC2(id2) = (VA(3,id2) + VA(4,id2)) - VC0(id2);
VC3(id3) = (VA(5,id3) + VA(6,id3)) - VC0(id3);


% Compute face gradients
dVdxE1 =  0.5.*( (VC1-VC0).*(yv(2,:)-yv(1,:)) + (VA(1,:)-VA(2,:)).*(yc1 - yc0) )./A1;
dVdyE1 = -0.5.*( (VC1-VC0).*(xv(2,:)-xv(1,:)) + (VA(1,:)-VA(2,:)).*(xc1 - xc0) )./A1;
dVdxE2 =  0.5.*( (VC2-VC0).*(yv(3,:)-yv(2,:)) + (VA(3,:)-VA(4,:)).*(yc2 - yc0) )./A2;
dVdyE2 = -0.5.*( (VC2-VC0).*(xv(3,:)-xv(2,:)) + (VA(3,:)-VA(4,:)).*(xc2 - xc0) )./A2;
dVdxE3 =  0.5.*( (VC3-VC0).*(yv(1,:)-yv(3,:)) + (VA(5,:)-VA(6,:)).*(yc3 - yc0) )./A3;
dVdyE3 = -0.5.*( (VC3-VC0).*(xv(1,:)-xv(3,:)) + (VA(5,:)-VA(6,:)).*(xc3 - xc0) )./A3;

dVdxC0 = (A1.*dVdxE1 + A2.*dVdxE2 + A3.*dVdxE3)./(A1+A2+A3);
dVdyC0 = (A1.*dVdyE1 + A2.*dVdyE2 + A3.*dVdyE3)./(A1+A2+A3);

dVdxC1 = dVdxC0(E1); dVdxC2 = dVdxC0(E2); dVdxC3 = dVdxC0(E3);
dVdyC1 = dVdyC0(E1); dVdyC2 = dVdyC0(E2); dVdyC3 = dVdyC0(E3);

% Build weights used in limiting
g1 = (dVdxC1.^2 + dVdyC1.^2); g2 = (dVdxC2.^2 + dVdyC2.^2); g3 = (dVdxC3.^2 + dVdyC3.^2);

epse = 1e-10; fac = g1.^2 + g2.^2 + g3.^2;
w1 = (g2.*g3+epse)./(fac+3*epse);
w2 = (g1.*g3+epse)./(fac+3*epse);
w3 = (g1.*g2+epse)./(fac+3*epse);

% Limit gradients
LdVdxC0 = w1.*dVdxC1 + w2.*dVdxC2 + w3.*dVdxC3;
LdVdyC0 = w1.*dVdyC1 + w2.*dVdyC2 + w3.*dVdyC3;

% the limiting result
hlim = dx.*(ones(mesh.Shape.nNode,1)*LdVdxC0) + ...
    dy.*(ones(mesh.Shape.nNode,1)*LdVdyC0) + ...
    ones(mesh.Shape.nNode,1)*VC0;
end% func