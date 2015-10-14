function error = LinConError(mesh, c)
% 1D Linear Convection Test Case
% Purpose: calculate error, 0.6s
% OUTPUT:
%   error   - [L1, L2, Linf]
Cout = zeros(size(mesh.x));
flag = (mesh.x < 0.7) & (mesh.x > 0.6);
% temp_c = c(flag);
temp_x = mesh.x(flag);

temp_c = sin(10*pi*(temp_x-0.6));
Cout(flag) = temp_c;

error(1) = sum(sum(abs(c - Cout)))./ mesh.nNode;
error(2) = sqrt( sum(sum( (c-Cout).^2 ))./mesh.nNode );
error(3) = max( max( abs(c - Cout) ) );
end% func