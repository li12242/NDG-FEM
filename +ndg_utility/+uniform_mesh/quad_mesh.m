function [K,EToV,Nv,VX,VY,EToBS,EToR] = quad_mesh(Mx, My, ...
    xmin, xmax, ymin, ymax, bc_type)
%QUAD_MESH Generation uniform rectangle mesh, where the domain is
%   of [xmin, xmax]x[ymin, ymax].
% Parameters:
%   Nx     - No of points along x
%   Ny     - No of points along y
%   xmin   - start position of x
%   xmax   - end position of x
%   ymin   - start position of y
%   ymax   - end position of y
% Returns:
%   EToV   - element to vertex list
%   VX     - vertex coordinate
%   VY     - vertex coordinate
% 
% Usages
%   [E, X, Y] = MeshGenRectangle2D(2, 3, -1.0, 1.0, -1.0, 1.0)
% 

%% Parameters
Nx = Mx + 1; % number of elements along x coordinate
Ny = My + 1;
K = Mx * My;
Nv = Nx * Ny;
EToR = int8( ones(K, 1) )*ndg_lib.mesh_type.Normal;

%% Define vectex
% The vertex is sorted along x coordinate. (x coordinate counts first)
VX   = linspace(xmin, xmax, Nx) ;
VY   = linspace(ymin, ymax, Ny)';
VX   = repmat(VX, 1, Ny) ;
VY   = repmat(VY, 1, Nx)';
VX   = VX(:);
VY   = VY(:);

%% Define EToV
% The element is conuting along x coordinate
EToV = zeros(4, Mx*My);
ind = 1;
for i = 1:My
    for j = 1:Mx
        % vertex index
        v1 = Nx*(i-1) + j;
        v2 = Nx*(i-1) + j + 1;
        v3 = Nx*i + j;
        v4 = Nx*i + j + 1;
        % Counterclockwise
        EToV(:, ind)=[v1, v2, v4, v3]';
        ind = ind +1;
    end% for
end% for

[ EToBS ] = ndg_utility.uniform_mesh.uniform_bc( Mx, My, EToV, bc_type );

end% func