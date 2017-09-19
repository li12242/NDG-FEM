function [ dr, ds, dt ] = derivative_orthogonal_func(obj, N, ind, r, s, t)
%DERIVATIVE_ORTHOGONAL_FUNC Summary of this function goes here
%   Detailed explanation goes here

[i, j, k] = trans_ind(N, ind);

[ f_quad ] = quad_orthgonal_func( i,j,r,s );
[ f_line ] = line_orthgonal_func( k,t );
[ dfdt ] = line_deri_orthgonal_func( k,t );
[ dfdr, dfds ] = quad_deri_orthgonal_func( i,j,r,s );

dr = dfdr.*f_line;
ds = dfds.*f_line;
dt = f_quad.*dfdt;
end

