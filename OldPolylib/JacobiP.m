function [p] = JacobiP(r,alpha,beta,n)
% normalized Jacobi polynamial
hpath = fullfile('/Users/mac/Documents/Model/NDG-FEM/+Polylib/libpolylib.h');
libpath = fullfile('/Users/mac/Documents/Model/NDG-FEM/+Polylib/libpolylib.dyld');
loadlibrary(libpath,hpath);

np = int32(numel(r));
n = int32(n);
% alpha = 0; beta =0;
p = zeros(1,np);


[~, p] = calllib('libpolylib', 'jacobiP', np, r, p, n, alpha, beta);
unloadlibrary libpolylib
end