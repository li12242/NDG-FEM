function [A] = PoissonCstab1D(Element, Mesh)

% function [A] = PoissonCstab1D();
% Purpose: Set up symmetric Poisson matrix with estabilized central fluxes

% Globals1D;
A = zeros(Mesh.K*Element.Np); g = zeros(Mesh.K*Element.Np,1);

% Build matrix -- one column at a time
for i=1:Mesh.K*Element.Np
    g(i) = 1.0;
    gmat = reshape(g,Element.Np,Mesh.K);
    Avec = PoissonCstabRHS1D(gmat, Element, Mesh);
    A(:,i) = reshape(Avec,Mesh.K*Element.Np,1);
    g(i)=0.0;
end
return
