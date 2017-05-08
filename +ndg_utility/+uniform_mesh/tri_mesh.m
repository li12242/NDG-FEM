function [K,EToV,Nv,VX,VY,EToBS,EToR] = tri_mesh(Nx, Ny, xmin, xmax, ymin, ymax, flag)
%TRI_MESH Summary of this function goes here
%   Generation uniform triangle mesh, 
%   the domain is of [xmin, xmax]x[ymin, ymax].
%
% INPUT
%    Nx     No. of points along x coordinate
%    Ny     No. of points along x coordinate
%    xmin   min x coordinate
%    xmax   max x coordinate
%    ymin   min y coordinate
%    ymax   max y coordinate
%    flag   = partination type [flag=false, "\", flag=true, "/"]
%
% OUTPUT
%    X = x coordinate
%    Y = y coordinate
%    EToV = index of vertices in element
%
% EXAMPLE USAGE
%    [EToV,X,Y] = MeshGenTriangle2D(2, 2, -1.0, 1.0, -1.0, 1.0, true)
%
% Author(s)
%    li12242 Tianjin University

%% Parameters
Mx   = Nx - 1; % number of elements along x coordinate
My   = Ny - 1;

%% Define vertex
% The vertex is sorted along x coordinate. (x coordinate counts first)
VX   = linspace(xmin, xmax, Nx) ;
VY   = linspace(ymin, ymax, Ny)'; 
VX   = repmat(VX, 1, Ny) ;
VY   = repmat(VY, 1, Nx)'; 
VX   = VX(:);
VY   = VY(:);

%% Define EToV
% The element is conuting along x coordinate
EToV = zeros(2*Mx*My,3);
for i = 1:My % each row
    for j = 1:Mx
        % element index
        ind1 = 2*Mx*(i-1)   +j;
        ind2 = 2*Mx*(i-1)+Mx+j;
        % vertex index
        v1   = Nx*(i-1) + j;
        v2   = Nx*(i-1) + j + 1;
        v3   = Nx*i + j;
        v4   = Nx*i + j + 1;
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