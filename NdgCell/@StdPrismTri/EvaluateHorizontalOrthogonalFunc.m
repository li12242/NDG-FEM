function [ f ] = EvaluateHorizontalOrthogonalFunc( obj, N1, td, r, s )
%> @brief Value of the horizontal (triangle) orthgonal basis.

td1 = mod( td - 1, obj.Nph ) + 1;
% project the coordinate (r,s) in triangle to (a,b) in square
[ a, b ] = ndgcell.rstoab( r, s );
% linear index to two indexes (i,j) on the collapsed square
[ i, j ] = ndgcell.trans_ind( N1, td1 );
[ f ] = ndgcell.simplex2DP( a, b, i, j );
end
