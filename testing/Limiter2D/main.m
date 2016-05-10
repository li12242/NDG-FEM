function main

mesh = quadSolver(N, M);
% mesh = triSolver(N, M);
var = ConvectionInit(mesh);

end% func


function var = ConvectionInit(mesh)
% initial condition

% left up, smooth part
sigma = 125*1e3/33^2; 
xc = -1/4; yc = 1/4;
flag = (mesh.x<=0) & (mesh.y>0);
var(flag) = exp(-sigma.*( (mesh.x(flag) - xc).^2 + (mesh.y(flag) - yc).^2) );

% left right, discontinuous part

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