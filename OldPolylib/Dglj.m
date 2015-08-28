function [D] = Dglj(r)
% Compute the Derivative Matrix and its transpose associated
%     with the Gauss-Lobatto-Jacobi zeros.
% 
% Compute the derivative matrix, associated with the n_th order 
%     Lagrangian interpolants through the np Gauss-Lobatto-Jacobi 
%     points z such that
%     $\frac{du}{dz}(z[i]) =  \sum_{j=0}^{np-1} D[i*np+j] u(z[j])$
%     $D[i*np+j] = \frac{\partial l_j}{\partial z} \right|_{z=z_i}$
% 
% Input:
%       double[np]  r
% Output:
%       double[np x np] D
% Usages: D=Polylib.Dglj(r)

hpath = fullfile('/Users/mac/Documents/Model/NDG-FEM/+Polylib/libpolylib.h');
libpath = fullfile('/Users/mac/Documents/Model/NDG-FEM/+Polylib/libpolylib.dyld');
loadlibrary(libpath,hpath);

[row,~] = size(r);
if (row ~= 1)
    r = r';
end
assert(size(r,1)==1, 'input should be a vector');

np = int32(numel(r));
D = zeros(1, np*np);

[D,~] = calllib('libpolylib', 'Dglj', D, D, r, np, 0, 0);
D = reshape(D, np, np);
unloadlibrary libpolylib
end