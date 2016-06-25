function A = SteadyConvectionCoeffMatrix(mesh)
% set up symmetric matrix

A = zeros(mesh.nNode, mesh.nNode);
g = zeros(size(mesh.x));

% Build matrix -- one column at a time
for i = 1:mesh.nNode
    g(i) = 1;
    
    Avec = SteadyConvectionRHS(mesh, g);
    A(:, i) = Avec(:);
    g(i) = 0;
end% for

end% func