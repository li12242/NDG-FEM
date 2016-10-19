% test FindLocCell
N   = 2;
tri  = StdRegions.Triangle(N);
quad = StdRegions.Quad(N);

% mesh
Ne = 2;
Nx = Ne+1; 
Ny = Ne+1;
rmin = -1;
rmax =  1;

[EToV,VX,VY] = Utilities.Mesh.MeshGenTriangle2D(Nx,Ny,rmin,rmax,rmin,rmax,1);
triMesh      = MultiRegions.RegionTri(tri, EToV, VX, VY);
triVerList   = tri.getVertexNodeList;

[EToV,VX,VY] = Utilities.Mesh.MeshGenRectangle2D(Nx,Ny,rmin,rmax,rmin,rmax);
quadMesh     = MultiRegions.RegionQuad(quad, EToV, VX, VY);
quadVerList  = quad.getVertexNodeList;

% error tolerance
tol = 1e-12;

%% Test 1: triangular mesh
xp = [-0.5,  0.25, 0.25, 0.75];
yp = [0.75, -0.25, 0.75, -0.75];
eleInd = [5, 2, 6, 4];
flag   = Utilities.PostProcess.FindLocCell_Mex...
    (triMesh.x,triMesh.y,triVerList,xp,yp);
            
assert( abs(flag(1) - eleInd(1)) <= tol)
assert( abs(flag(2) - eleInd(2)) <= tol)
assert( abs(flag(3) - eleInd(3)) <= tol)
assert( abs(flag(4) - eleInd(4)) <= tol)

%% Test 2: quadrilateral mesh
xp = [-0.5,  0.25, 0.25, 0.75];
yp = [0.75, -0.25, 0.75, -0.75];
eleInd = [3, 2, 4, 2];
flag   = Utilities.PostProcess.FindLocCell_Mex...
    (quadMesh.x,quadMesh.y,quadVerList,xp,yp);
            
assert( abs(flag(1) - eleInd(1)) <= tol)
assert( abs(flag(2) - eleInd(2)) <= tol)
assert( abs(flag(3) - eleInd(3)) <= tol)
assert( abs(flag(4) - eleInd(4)) <= tol)

%% Test 3: boundary nodes
xp = [0, 0, 1, 0];
yp = [-0.25, .25, 0.25, .25];
eleInd = [1, 3, 4, 3];
flag   = Utilities.PostProcess.FindLocCell_Mex...
    (quadMesh.x,quadMesh.y,quadVerList,xp,yp);
            
assert( abs(flag(1) - eleInd(1)) <= tol)
assert( abs(flag(2) - eleInd(2)) <= tol)
assert( abs(flag(3) - eleInd(3)) <= tol)