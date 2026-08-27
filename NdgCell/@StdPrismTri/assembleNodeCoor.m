function assembleNodeCoor( obj, Nh, Nz )
%> @brief Assemble the interpolation node coordinates of the triangular prism,
%> as the tensor product of triangle nodes (ndgcell.triNodeCoor) in the
%> horizontal direction and LGL nodes in the vertical direction.

[ Nph, r1, s1 ] = ndgcell.triNodeCoor( Nh );
Npz = Nz + 1;
[ t1, ~ ] = zwglj( Npz );

Np = Nph * Npz;
r = repmat(r1, 1, Npz);
s = repmat(s1, 1, Npz);
t = repmat(t1', Nph, 1 );

obj.Np = Np;
obj.Nph = Nph;
obj.Npz = Npz;

obj.r1 = r1;
obj.s1 = s1;
obj.t1 = t1;

obj.r = r(:);
obj.s = s(:);
obj.t = t(:);
end
