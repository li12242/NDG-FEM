function [ f ] = orthogonal_func(obj, N, ind, r, s, t)
%ORTHOGONAL_FUNC Summary of this function goes here
%   Detailed explanation goes here

[i, j, k] = trans_ind(N, ind);

[ f_quad ] = quad_orthgonal_func( i, j, r, s );
[ f_line ] = line_orthgonal_func( k, t );
[ f ] = f_quad .* f_line;
end

