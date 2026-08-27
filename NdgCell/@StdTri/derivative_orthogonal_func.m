function [ dr, ds, dt ] = derivative_orthogonal_func( obj, N, ind, r, s, t )
%> @brief Derivatives of the orthgonal basis on the reference triangle.

% project to square.
[ a, b ] = ndgcell.rstoab( r, s );
% transform the index to two indexes.
[ i, j ] = ndgcell.trans_ind( N, ind );
% calculate the derivative function values.
[ dr, ds ] = ndgcell.deriSimplex2DP( a, b, i, j );
dt = zeros(size(dr));
end
