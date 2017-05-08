function [K,EToV,Nv,VX,VY,EToBS,EToR] = quad_mesh(Nx, Ny, xmin, xmax, ymin, ymax)
%NDG_CELL Summary of this function goes here
%   Generation uniform rectangle mesh, 
%   The domain is of [xmin, xmax]x[ymin, ymax].
% Input:
%   Nx     - No of points along x
%   Ny     - No of points along y
%   xmin   - start position of x
%   xmax   - end position of x
%   ymin   - start position of y
%   ymax   - end position of y
% Output:
%   EToV   - element to vertex list
%   VX     - vertex coordinate
%   VY     - vertex coordinate
% 
% Usages
%   [E, X, Y] = MeshGenRectangle2D(2, 3, -1.0, 1.0, -1.0, 1.0)
% 

% Parameters
Mx   = Nx - 1; % number of elements along x coordinate
My   = Ny - 1;

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
EToV = zeros(Mx*My, 4);
ind = 1;
for i = 1:My
    for j = 1:Mx
        % vertex index
        v1   = Nx*(i-1) + j;
        v2   = Nx*(i-1) + j + 1;
        v3   = Nx*i     + j;
        v4   = Nx*i     + j + 1;
        % Counterclockwise
        EToV(ind,:)=[v1, v2, v4, v3];
        ind = ind +1;
    end% for
end% for

end% func