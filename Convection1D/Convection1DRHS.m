function rhsu = Convection1DRHS(mesh, u, a)

line = mesh.Shape;

% form field differences at faces
alpha=0;
% du = zeros(size(mesh.nx));
du = (u(mesh.vmapM)-u(mesh.vmapP)).*...
    (a.*mesh.nx-(1-alpha)*abs(a.*mesh.nx))/2;

% impose boundary condition


% compute right hand sides of the semi-discrete PDE
rhsu = -mesh.rx.*(line.Dr*(a.*u)) + line.invM*line.Mes*( du.*mesh.fScale );

end% func
