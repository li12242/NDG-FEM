function Mesh1D = BuildMaps1D(Mesh1D)

% function [vmapM, vmapP, vmapB, mapB] = BuildMaps1D
% Purpose: Connectivity and boundary tables for nodes given in the K # of elements,
% 	       each with N+1 degrees of freedom.

% Globals1D;
NODETOL = 10e-5;

K = Mesh1D.K;  Nfp = Mesh1D.Element.Nfp;
Np = Mesh1D.Element.Np; Nfaces = Mesh1D.Element.Nfaces;

% number volume nodes consecutively
nodeids = reshape(1:K*Np, Np, K);
vmapM   = zeros(Nfp, Nfaces, K);
vmapP   = zeros(Nfp, Nfaces, K);

for k1=1:K
  for f1=1:Nfaces
    % find index of face nodes with respect to volume node ordering
    vmapM(:,f1,k1) = nodeids(Mesh1D.Element.Fmask(:,f1), k1);
  end
end

for k1=1:K
  for f1=1:Nfaces
    % find neighbor
    k2 = Mesh1D.EToE(k1,f1); f2 = Mesh1D.EToF(k1,f1);
    
    % find volume node numbers of left and right nodes 
    vidM = vmapM(:,f1,k1); vidP = vmapM(:,f2,k2);
    
    x1  = Mesh1D.x(vidM); x2  = Mesh1D.x(vidP);
    
    % Compute distance matrix
    D = (x1 -x2 ).^2;
    if (D<NODETOL) vmapP(:,f1,k1) = vidP; end;
  end
end

Mesh1D.vmapP = vmapP(:); Mesh1D.vmapM = vmapM(:);

% Create list of boundary nodes
Mesh1D.mapB = find(vmapP==vmapM); Mesh1D.vmapB = vmapM(Mesh1D.mapB);

% Create specific left (inflow) and right (outflow) maps
Mesh1D.mapI = 1; Mesh1D.mapO = K*Nfaces; Mesh1D.vmapI = 1; Mesh1D.vmapO = K*Np;
return
