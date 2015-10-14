% Driver script for solving the 1D Poisson equation
% Globals1D;

% Polynomial order used for approximation 
N =4;

Element = Element1D(N);

% Read in Mesh
[Nv, VX, K, EToV] = MeshGen1D(0,2*pi,10);

Mesh = Mesh1D(Element, EToV, VX);

% Initialize solver and construct grid and metric
% StartUp1D;

% Set RHS
f = -Mesh.J.*((Element.invV'*Element.invV)*sin(Mesh.x));

% Set up operator
[A] = PoissonCstab1D(Element, Mesh);

% Solve Problem
solvec = A\f(:);
u = reshape(solvec,Element.Np,Mesh.K);
