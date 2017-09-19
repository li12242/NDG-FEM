function [ fval ] = quad_orthgonal_func( i, j, r, s )
%QUAD_ORTHGONAL_FUNC Summary of this function goes here
%   Detailed explanation goes here

fval = Polylib.JacobiP(r(:), 0, 0, i).*Polylib.JacobiP(s(:), 0, 0, j);
end

