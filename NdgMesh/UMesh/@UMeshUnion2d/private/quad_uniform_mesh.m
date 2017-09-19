function [Nv,VX,VY,K,EToV,EToR,EToBS] ...
    = quad_uniform_mesh( xlim, ylim, Mx, My, bc_type )
%QUAD_UNIFORM_MESH Generate informations for the uniform quadrilateral mesh
%
%   The mesh domain is a rectangle controlled by 'xlim' and 'ylim' for the 
%   length on x and y coordinates. The number of elements on each
%   coordinates is set by 'Mx' and 'My'. All the quadrilateral elements are 
%   rectangles and the number is (Mx*My). 
%
%   The nodes are numbered from the left bottom corner and the loop through
%   x coordinate is numbered first. The elements are also numbered from the
%   bottom to the top. 
% 
%   The input parameters includes:
%       xlim, ylim  - the ranges of the domain on x and y coordinates
%       Mx, My      - number of elements on the x and y coordinates
%       bctype      - The boundary condition of the sourth, north, west and
%                       east sides
%
%   The output parameters:
%       Nv          - number of the vertices
%       vx, vy      - coordinates of the vertices
%       K           - number of the elements
%       EToV        - the index of the vertex in each element (column)
%       EToR        - the element type id
%       EToBS       - edge type for each element (column)
%
% Author: li12242 Tianjin University
%

%% Parameters
Nx = Mx + 1; % number of elements along x coordinate
Ny = My + 1;
K = Mx * My;
Nv = Nx * Ny;
EToR = int8( ones(K, 1) )*ndg_lib.mesh_type.Normal;

%% Define vectex
% The vertex is sorted along x coordinate. (x coordinate counts first)
xmin = min(xlim); xmax = max(xlim);
ymin = min(ylim); ymax = max(ylim);
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

end

