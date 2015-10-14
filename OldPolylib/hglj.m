function [h] = hglj(i,z,zglj,np)
% Compute the value of the i th Lagrangian interpolant through the
% np Gauss-Lobatto-Jacobi points zgrj at the arbitrary location z.
% Input:
%       int     i
%       double  z
%       double[np]  zglj
%       int     np
% Output:
%       double  h
% 
% Usages: [h] = hglj(i,z,zglj,np)


hpath = fullfile('/Users/mac/Documents/Model/NDG-FEM/+Polylib/libpolylib.h');
libpath = fullfile('/Users/mac/Documents/Model/NDG-FEM/+Polylib/libpolylib.dyld');
loadlibrary(libpath,hpath);

i = int32(i);
np = int32(np);
assert(numel(zglj) == np, 'the number of nodal base is not correct');

[h,~] = calllib('libpolylib', 'hglj', i, z, zglj, np, 0, 0);
unloadlibrary libpolylib
end