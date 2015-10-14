function rhs = LinConRHS(mesh, c, u, time)
% 1D Linear Convection Test Case
% Purpose  : Evaluate RHS flux in strong form

line = mesh.Shape;
cs = LinConNumFlux(mesh,u ,c);
dc = (c(mesh.vmapM).*mesh.nx - cs);
% alpha = 0;
% dc = (c(mesh.vmapM)-c(mesh.vmapP)).*(u*mesh.nx-(1-alpha)*abs(u*mesh.nx))/2;

dc = LinConBC(mesh,c,dc,u, time);

rhs = - mesh.rx.*(line.Dr*(u.*c)) + line.invM*line.FaceMassMatrixSmall*(mesh.fScale.*dc);
end% func