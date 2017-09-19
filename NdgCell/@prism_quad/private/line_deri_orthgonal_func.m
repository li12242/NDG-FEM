function [ df ] = line_deri_orthgonal_func( ind, r )
%LINE_DERI_ORTHGONAL_FUNC Summary of this function goes here
%   Detailed explanation goes here

df = Polylib.GradJacobiP(r, 0, 0, ind-1);        
end

