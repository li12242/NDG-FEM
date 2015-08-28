function [dp] = GradJacobiP(r, alpha, beta, n)
hpath = fullfile('/Users/mac/Documents/Model/NDG-FEM/+Polylib/libpolylib.h');
libpath = fullfile('/Users/mac/Documents/Model/NDG-FEM/+Polylib/libpolylib.dyld');
loadlibrary(libpath,hpath);

np = int32(numel(r));
n = int32(n);
% alpha = 0; beta =0;
dp = zeros(1,np);

[~, dp] = calllib('libpolylib', 'GradjacobiP', np, r, dp, n, alpha, beta);
unloadlibrary libpolylib
end