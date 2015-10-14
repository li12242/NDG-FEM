function error = CompareDerivaMatrix(N)
%% Polylib
nModes = N+1;
[nodalPoints,~] = Polylib.zwglj(nModes);
D = Polylib.Dglj(nodalPoints);
%% NDGFEM
r = nodalPoints;
V  = Vandermonde1D(N, r);
D1 = Dmatrix1D(N, r, V);
%% compare
error = abs(D'-D1);

end