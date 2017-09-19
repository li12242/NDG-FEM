function [ P ] = tri_orthgonal_func( i,j,r,s )
% Evaluate 2D orthonormal polynomial on simplex at (a,b) of order (i,j).
% 

% project the coordinate from the triangle to the square
[ a,b ] = rstoab( r,s );

h1 = Polylib.JacobiP(a,0,0,i); 
h2 = Polylib.JacobiP(b,2*i+1,0,j);
P = sqrt(2.0)*h1.*h2.*(1-b).^i;
end

