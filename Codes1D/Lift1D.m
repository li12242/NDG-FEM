function [LIFT] = Lift1D(Element1D)

% function [LIFT] = Lift1D
% Purpose  : Compute surface integral term in DG formulation

% Globals1D;
Emat = zeros(Element1D.Np, Element1D.Nfaces*Element1D.Nfp);

% Define Emat
Emat(1,1) = 1.0; Emat(Element1D.Np,2) = 1.0;

% inv(mass matrix)*\s_n (L_i,L_j)_{edge_n}
LIFT = Element1D.V*(Element1D.V'*Emat);
return
