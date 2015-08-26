function Cout = LinConInit(mesh)
% 1D Linear Convection Initial Condition
% REFERENCE:
% [1] Khan A A, Lai W. Modeling Shallow Water Flows Using the 
%     Discontinuous Galerkin Method[M]. CRC Press, 2014. 49-52

Cout = zeros(size(mesh.x));
flag = mesh.x < 0.1;
% temp_c = c(flag);
temp_x = mesh.x(flag);

temp_c = sin(10*pi*temp_x);
Cout(flag) = temp_c;
end% fun