function [dr, ds, dt] = derivative_orthogonal_func(obj, N, ind, r, s, t)
%DERIVATIVE_ORTHOGONAL_FUNC Summary of this function goes here
%   Detailed explanation goes here

% transform the 
[ i, j, k ] = trans_ind( N,ind );
% project the coordinate from the triangle to the square
[ f_tri ] = tri_orthgonal_func( i,j,r,s );
[ f_lin ] = line_orthgonal_func(k, t);
[ dfdt ] = line_deri_orthgonal_func(k, t);
[ dfdr, dfds ] = tri_deri_orthgonal_func( i,j,r,s );

dr = dfdr.*f_lin;
ds = dfds.*f_lin;
dt = f_tri.*dfdt;
end
