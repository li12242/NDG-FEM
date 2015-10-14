function [nx] = Normals1D(Mesh1D)

% function [nx] = Normals1D
% Purpose : Compute outward pointing normals at elements faces

% Globals1D;
Nfp = Mesh1D.Element.Nfp;
Nfaces = Mesh1D.Element.Nfaces;

nx = zeros(Nfp*Nfaces, Mesh1D.K); 

% Define outward normals
nx(1, :) = -1.0; nx(2, :) = 1.0;

return
