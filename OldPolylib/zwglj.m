function [z,w] = zwglj(np)
% Gauss-Lobatto-Jacobi zeros and weights with end point at $z\in[-1, 1]$
% Input:    
%       int     np
% Output:
%       double  z
%       double  w
% 
% Usages: [z,w] = zwglj(np)

hpath = fullfile('/Users/mac/Documents/Model/NDG-FEM/+Polylib/libpolylib.h');
libpath = fullfile('/Users/mac/Documents/Model/NDG-FEM/+Polylib/libpolylib.dyld');
loadlibrary(libpath,hpath);

np = int32(np);
z = zeros(1,np);
w = zeros(1,np);

[z, w] = calllib('libpolylib', 'zwglj', z,w,np,0,0);
unloadlibrary libpolylib
end% func