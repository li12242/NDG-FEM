function [ P ] = simplex2DP( a, b, i, j )
%> @brief Evaluate 2D orthonormal polynomial on simplex at (a,b) of order (i,j).
% ======================================================================
%> This function is part of the NDG-FEM software.
% ======================================================================
h1 = JacobiP(a, 0, 0, i);
h2 = JacobiP(b, 2*i+1, 0, j);
P = sqrt(2.0)*h1.*h2.*(1-b).^i;
end% func
