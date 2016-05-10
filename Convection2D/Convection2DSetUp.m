function [mesh, var] = Convection2DSetUp(N, M)
% 2D convection problem
% dc/dt + d(uc)/dx + d(vc)/dy = 0
% Input:
%   N - degree of polynomial
%   M - No. of elements on each edge
% 

mesh = quadSolver(N, M);
% mesh = triSolver(N, M);
var = ConvectionInit(mesh);

w = 5*pi/6;
u = -w.*mesh.y; % flow rate in domain, [u, v]
v = w.*mesh.x;

FinalTime = 0.6;
filename = ['Convection2D_', num2str(N),'_',num2str(M),'.nc'];
outfile = CreateNetcdfFile(mesh, filename);

var = Convection2DSolver(mesh, var, FinalTime, u, v, outfile);
% postprocess(mesh, var)
end% func

function var = ConvectionInit(mesh)
% var = ones(size(mesh.x));
% var = mesh.x;
% xc = mean(mesh.x);
% left = xc < 0.5; right = xc > 0.5;
% var = sin(pi*mesh.x);%.*sin(2*pi*mesh.y);
% var = zeros(size(mesh.x));
% var(:,left) = 1; var(:,right) = 0;

sigma = 125*1e3/33^2; 
xc = 0; yc = 3/5;
var = exp(-sigma.*( (mesh.x - xc).^2 + (mesh.y - yc).^2) );
end% func

function mesh = triSolver(N, M)

% read triangle mesh
% [EToV, VX, VY, EToR, BC] = Utilities.Mesh.MeshReaderTriangle('Convection2D/mesh/triangle');

np = M+1;
[VX,VY,EToV] = Lmesh(np, -1, 1, 1);

tri = StdRegions.Triangle(N);
mesh = MultiRegions.RegionTri(tri, EToV, VX, VY);
end% func

function mesh = quadSolver(N, M)

% read triangle mesh
% [EToV, VX, VY, EToR, BC] = Utilities.Mesh.MeshReaderQuad('Convection2D/mesh/quad');

% uniform mesh
[EToV, VX, VY] = uniformRetangle(M+1);

temp = EToV(:, 3); 
EToV(:, 3) = EToV(:, 4);
EToV(:, 4) = temp;

quad = StdRegions.Quad(N);
mesh = MultiRegions.RegionQuad(quad, EToV, VX, VY);

end% func

function [X,Y,E] = Lmesh(n,start_coor,end_coor,flag)
% Mesh generater with triangle element
%
% DESCRIPTION
%   domain [start_coor, end_coor]x[start_coor, end_coor];
%   n - No. of element on each edge
%
%   index
%       x coordinate 
% 
%       1   2   3   4   5   ....            n  n+1
%       -----------------------------------------
%       |   |   |   |   |   |   |   |   |   |   |
% (n+2) ----------------------------------------- 2(n+1)
%
% INPUT
%    n      = No. of elements on each edge
%    start_coor = start coordinate
%    end_coor   = end coordinate
%    flag   = partination type [flag=0, "\", flag=1, "/"]
%
% OUTPUT
%    X = x coordinate
%    Y = y coordinate
%    E = index of vertices in element
%
% EXAMPLE USAGE
%    [X,Y,E] = Lmesh(21, 0, 1, 1)
%
% Author(s)
%    li12242 Tianjin University
%
% Reversion
%    v1.0 2014-12-8
%========================================================================== 

% info of point
x = linspace(start_coor, end_coor, n+1); y = linspace(end_coor, start_coor, n+1);
x = repmat(x, n+1, 1); y = repmat(y, n+1, 1);
x = x';
X = x(:); Y = y(:);

E = zeros(2*n^2,3);
for i = 1:n % each row
    upLayerNum = (n+1)*(i-1)+1:(n+1)*(i-1)+n+1;
    downLayerNum = upLayerNum + (n+1);

    if flag     %
        E(2*n*(i-1)+1:2*n*(i-1)+n,:) = [downLayerNum(1:n)', upLayerNum(2:n+1)', upLayerNum(1:n)'];
        E(2*n*(i-1)+n+1:2*n*i,:) = [downLayerNum(1:n)',downLayerNum(2:n+1)',upLayerNum(2:n+1)'];
        flag=~flag;
    else
        E(2*n*(i-1)+1:2*n*(i-1)+n,:) = [downLayerNum(2:n+1)', upLayerNum(2:n+1)', upLayerNum(1:n)'];
        E(2*n*(i-1)+n+1:2*n*i,:) = [downLayerNum(1:n)',downLayerNum(2:n+1)',upLayerNum(1:n)'];
        flag=~flag;
    end

end
end

function [EToV, VX, VY] = uniformRetangle(np)
m = np - 1; % elements on each row

VX = linspace(-1, 1, np); VX = repmat(VX, 1, np);
% dx = [-1, 0, -1, 0, -1];
% for i = 1:np
%     VX((i-1)*np+1: i*np) = VX((i-1)*np+1: i*np) + dx(i);
% end
VY = linspace(1, -1, np)'; VY = repmat(VY, 1, np)'; VY = VY(:);

EToV = zeros(m^2, 4);
for irow = 1:m
    for icol = 1:m
        index = (irow-1)*m + icol;
        EToV(index, :) = [np*(irow)+icol, np*(irow)+icol+1, ...
            np*(irow-1)+icol, np*(irow-1)+icol+1];
%         EToV(index, :) = [np*(irow-1)+icol, np*(irow-1)+icol+1,...
%             np*(irow)+icol, np*(irow)+icol+1];
    end
end

end% func

function postprocess(mesh, var)
plot3(mesh.x(mesh.vmapP),mesh.y(mesh.vmapP), var(mesh.vmapM))
end% func

