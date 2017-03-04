function rhs = SteadyConvectionRHS(mesh, u)
% right hands of equation
% 

us = zeros( size(u(mesh.vmapM)) );
us(1, :) = u(mesh.vmapP(1, :));
us(2, :) = u(mesh.vmapM(2, :));

us = us.*mesh.nx;

rhs = (mesh.Shape.Mes * us) - mesh.J.*( mesh.rx .*( mesh.Shape.Dr'*(mesh.Shape.M*u) ));
end% func