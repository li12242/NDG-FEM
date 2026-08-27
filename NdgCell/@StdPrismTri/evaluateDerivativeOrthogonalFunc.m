function [dr, ds, dt] = evaluateDerivativeOrthogonalFunc( obj, N1, N2, td, r, s, t)
%> @brief Derivatives of the prism orthgonal basis, as products of the
%> horizontal (triangle) and vertical (line) modes.

td1 = mod( td-1, obj.Nph ) + 1;
td2 = ceil( td / obj.Nph );

% horizontal part: triangle mode (td1) and its derivatives
[ a, b ] = ndgcell.rstoab( r, s );
[ i, j ] = ndgcell.transInd( N1, td1 );
[ drh, dsh ] = ndgcell.deriSimplex2DP( a, b, i, j );
% vertical mode value
[ ft ] = JacobiP(t, 0, 0, td2-1);
% multiply the vertical polynomial
dr = drh .* ft;
ds = dsh .* ft;

% vertical part: derivative of the line mode times the horizontal mode
[ dt ] = GradJacobiP(t, 0, 0, td2-1);
[ frs ] = ndgcell.simplex2DP( a, b, i, j );
dt = dt .* frs;
end
