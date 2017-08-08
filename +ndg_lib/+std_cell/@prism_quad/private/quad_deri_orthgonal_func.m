function [ fdr, fds ] = quad_deri_orthgonal_func( i, j, r, s )
%QUAD_DERI_ORTHGONAL_FUNC Summary of this function goes here
%   Detailed explanation goes here

fdr = Polylib.GradJacobiP(r(:), 0, 0, i).*Polylib.JacobiP(s(:), 0, 0, j);
fds = Polylib.JacobiP(r(:), 0, 0, i).*Polylib.GradJacobiP(s(:), 0, 0, j);
end

