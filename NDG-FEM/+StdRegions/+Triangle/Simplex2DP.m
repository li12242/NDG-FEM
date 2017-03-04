function [P] = Simplex2DP(a,b,i,j)
% Evaluate 2D orthonormal polynomial on simplex at (a,b) of order (i,j).
% 
h1 = Polylib.JacobiP(a,0,0,i); h2 = Polylib.JacobiP(b,2*i+1,0,j);
P = sqrt(2.0)*h1.*h2.*(1-b).^i;
end