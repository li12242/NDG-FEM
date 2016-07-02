function [Fhs, Fqx, Fqy] = SWENumFlux2d(phys, mesh, h, qx, qy, dryEleFlag)
% Calculate the numerical flux
% Input:
%   phys - strucutre variable, it contains
%       |
%       | - 
% 
%   mesh - mesh object
%   h    - water depth
%   qx   - water flux on x coordinate
%   qy   - water flux on y coordinate
%   dryEleFlag - flag for dry element 
%
% Output:
%   [Fhs, Fqxs, Fqys] - numerical flux of each eqs
% 
%% Usages
% 
%   shape        = StdRegions.Triangle(1);
%   [VX,VY,EToV] = Utilities.Mesh.MeshGenTriangle2D(1,-1,1,0);
%   mesh         = MultiRegions.RegionTri(shape, EToV, VX, VY);
%   h   = [0, 1e-4, 1e-3; 1, 1e-4, 1]';
%   qx  = [0, 1e-3, 1e-3; 0, 1e-3, 1e-4]';
%   qy  = [0, 1e-3, 2e-3; 1, 1e-3, 1e-3]';
%   dryEleFlag = [true, false];
% 
%   phys.gra = 9.81;
%   phys.mesh = mesh;
%   [Fhs, Fqxs, Fqys] = SWENumFlux2d(phys, mesh, h, qx, qy, dryEleFlag)
% 

%% Get wet/dry status of each nodes
% extend the elements' wet/dry status to each nodes
% dryNodeFlag = repmat(dryEleFlag, mesh.Shape.nNode, 1);
% minDepth    = phys.minDepth;
% dryNodeFlag(h <= minDepth) = 1;
% dryM        = dryNodeFlag(mesh.vmapM);
% dryP        = dryNodeFlag(mesh.vmapP);

%% Rotational invariance
% The numerical flux is calculated by the rotational invariance property of
% SWE. The rotational and its inverse matrix is 
% 
% $T = \left[ \begin{array}{ccc} 1&0& 0 \cr 0&n_x&n_y \cr 0&-n_y&n_x 
% \end{array} \right] \quad T^{-1} = \left[ \begin{array}{ccc} 1&0&0 
% \cr 0&n_x&-n_y \cr 0&n_y&n_x \end{array} \right]$
% 
% and the numerical flux is calculated by 
% 
% $\bf{F^*} = F^* \cdot n_x + G^* \cdot n_y = T^{-1}\cdot F^{LF}(T\cdot U)$
% 
% while the $F^{LF}$ is the Lax-Friedrichs numerical flux function. It is
% able to use other numerical function, e.g. HLL and HLLC.
% Refer to Lai and Khan (2012) for more details.
% 


%% Rotate the primary variable
% 
% $Q=T\cdot U=\left(\begin{array}{c}h\cr q_x\cdot n_x+q_y\cdot n_y \cr
% q_x\cdot -n_y +q_y\cdot n_x\end{array}\right)$
% 

QxM =  qx(mesh.vmapM).*mesh.nx + qy(mesh.vmapM).*mesh.ny;
QyM = -qx(mesh.vmapM).*mesh.ny + qy(mesh.vmapM).*mesh.nx;

QxP =  qx(mesh.vmapP).*mesh.nx + qy(mesh.vmapP).*mesh.ny;
QyP = -qx(mesh.vmapP).*mesh.ny + qy(mesh.vmapP).*mesh.nx;

hM  = h(mesh.vmapM);
hP  = h(mesh.vmapP);

% Lax-Friedrichs flux function
% [Fhs, Fqxs, Fqys] = LLFFlux(phys, hM, hP, QxM, QxP, QyM, QyP, dryM, dryP);
% HLL flux function
[Fhs, Fqxs, Fqys] = HLLFlux(phys, hM, hP, QxM, QxP, QyM, QyP, dryEleFlag);

% Rotate invariance
Fqx = Fqxs.*mesh.nx  - Fqys.*mesh.ny;
Fqy = Fqxs.*mesh.ny  + Fqys.*mesh.nx;

end% func