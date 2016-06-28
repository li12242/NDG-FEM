%% MeshGenRectangle2D
% Generation uniform mesh of rectangle elements.
% The domain is of [xmin, xmax]x[ymin, ymax], 

function [EToV, VX, VY] = MeshGenRectangle2D(nx, ny, xmin, xmax, ymin, ymax)
% Input:
%   nx     - No of vertice along x
%   ny     - No of vertice along y
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
mx   = nx - 1; % number of elements along x coordinate
my   = ny - 1;

%% Define vectex
% The vertex is sorted along x coordinate. (x coordinate counts first)
VX   = linspace(xmin, xmax, nx) ;
VY   = linspace(ymin, ymax, ny)'; 
VX   = repmat(VX, 1, ny) ;
VY   = repmat(VY, 1, nx)'; 
VX   = VX(:);
VY   = VY(:);

%% Define EToV
% The element is conuting along x coordinate
EToV = zeros(mx*my, 4);
for i = 1:my
    for j = 1:mx
        % element index
        ind  = (i-1)*mx + j;
        % vertex index
        v1   = nx*(i-1) + j;
        v2   = nx*(i-1) + j + 1;
        v3   = nx*i     + j;
        v4   = nx*i     + j + 1;
        % Counterclockwise
        EToV(ind,:)=[v1, v2, v4, v3];
    end% for
end% for

end% func