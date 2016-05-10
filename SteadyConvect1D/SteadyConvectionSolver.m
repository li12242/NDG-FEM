function u = SteadyConvectionSolver(mesh, f)
% set up and solve the equation system
% Input:
%   mesh - mesh object
%   f - source term
% Output:
%   u - unknown variable

f = mesh.J.*(mesh.Shape.M*f);

%% set up and solve global matrix coeffcient

% get system global matrix coefficient
A = SteadyConvectionCoeffMatrix(mesh);

% boundary condition
u0 = 1; M = 1e8;
A(1, 1) = M; f(1) = u0*M;

solvec = A\f(:);
u = reshape(solvec, size(mesh.x) );

end% func