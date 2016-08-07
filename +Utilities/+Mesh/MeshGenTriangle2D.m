%% MeshGenTriangle2D
% Generate uniform mesh with triangle elements
% 
function [VX,VY,EToV] = MeshGenTriangle2D(nx, ny, xmin, xmax, ymin, ymax, flag)
% DESCRIPTION
%   domain [start_coor, end_coor]x[start_coor, end_coor];
%   ne - No. of element on each edge
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
%    ne      = No. of elements on each edge
%    start_coor = start coordinate
%    end_coor   = end coordinate
%    flag   = partination type [flag=false, "\", flag=true, "/"]
%
% OUTPUT
%    X = x coordinate
%    Y = y coordinate
%    E = index of vertices in element
%
% EXAMPLE USAGE
%    [X,Y,E] = MeshGenTriangle2D(2, 2, -1.0, 1.0, -1.0, 1.0, true)
%
% Author(s)
%    li12242 Tianjin University

%% Parameters
mx   = nx - 1; % number of elements along x coordinate
my   = ny - 1;

%% Define vertex
% The vertex is sorted along x coordinate. (x coordinate counts first)
VX   = linspace(xmin, xmax, nx) ;
VY   = linspace(ymin, ymax, ny)'; 
VX   = repmat(VX, 1, ny) ;
VY   = repmat(VY, 1, nx)'; 
VX   = VX(:);
VY   = VY(:);

%% Define EToV
% The element is conuting along x coordinate
EToV = zeros(2*mx*my,3);
for i = 1:my % each row
    for j = 1:mx
        % element index
        ind1 = 2*mx*(i-1)   +j;
        ind2 = 2*mx*(i-1)+mx+j;
        % vertex index
        v1   = nx*(i-1) + j;
        v2   = nx*(i-1) + j + 1;
        v3   = nx*i     + j;
        v4   = nx*i     + j + 1;
        % Counterclockwise
        if flag % '/' divided
            EToV(ind1,:) = [v1, v4, v3];
            EToV(ind2,:) = [v1, v2, v4];
        else    % '\' divided
            EToV(ind1,:) = [v1, v2, v3];
            EToV(ind2,:) = [v2, v4, v3];
        end% if
    end
end
end