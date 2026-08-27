function [ fval ] = orthogonal_func( obj, N, td, r, s, t )
%> @brief Get the td-th orthgonal function value at the coordinate (r,s,t).

% project the coordinate (r,s) in triangle to (a,b) in square
[ a, b ] = ndgcell.rstoab( r, s );
% linear index to two indexes (i,j) on the collapsed square
[ i, j ] = ndgcell.trans_ind( N, td );
[ fval ] = ndgcell.simplex2DP( a, b, i, j );
end
