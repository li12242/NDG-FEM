% test for Riemann problem
tol = 1e-12;

%% Quad: N=1
% std element
N = 1;
shape = StdRegions.Quad(N);

% quadrilateral mesh
Nx = 5; 
Ny = 3;
rmin = -1;
rmax = 1;
smin = -1;
smax = 1;
[EToV, VX, VY] = Utilities.Mesh.MeshGenRectangle2D...
    (Nx,Ny,rmin,rmax,smin,smax);
mesh = MultiRegions.RegionQuad(shape, EToV, VX, VY);

% initialization
h = zeros(size(mesh.x));
ind = [1,5];
h(:, ind) = 10;
ind = [8,4];
h(:, ind) = 2;
ind = [6,2];
h(:, ind) = repmat([10.2007, 9.5985, 10.2007, 9.5985]', 1, 2);
ind = [7,3];
h(:, ind) = repmat([2.4015, 1.7993, 2.4015, 1.7993]', 1, 2);

hlim = Utilities.Limiter.Limiter2D.VB2d_VA(mesh, h);

%% Tri: N=1