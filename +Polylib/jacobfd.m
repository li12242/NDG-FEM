function [p,dp] = jacobfd(r, alpha, beta, n)
% Returns value and derivative of Jacobi poly. at point z
% 
% This function returns the vectors poly_in and poly_d
%     containing the value of the $ n^th $ order Jacobi polynomial
%     $ P^{\alpha,\beta}_n(z) \alpha > -1, \beta > -1 $ and its
%     derivative at the np points in z[i]
% Input:
%       int         n   nth order Jacobi polynomial
%       double[np]  r   np points $\in$ [-1,1]
% Output:
%       
% 

hpath = fullfile('/Users/mac/Documents/Model/NDG-FEM/+Polylib/libpolylib.h');
libpath = fullfile('/Users/mac/Documents/Model/NDG-FEM/+Polylib/libpolylib.dyld');
loadlibrary(libpath,hpath);

np = int32(numel(r));
n = int32(n);
% alpha = 0; beta =0;
p = zeros(1,np); dp = zeros(1,np);


[~, p, dp] = calllib('libpolylib', 'jacobfd', np, r, p, dp, n, alpha, beta);
unloadlibrary libpolylib

end