function rhsu = Convection1d_RHS(phys, var)

mesh = phys.mesh;
line = mesh.Shape;
a = phys.u;

mapI = 1;
vmapI = 1;

dT = ( phys.var(mesh.vmapM)-phys.var(mesh.vmapP) );
Tin = 0; dT(mapI) = 2.0*(var(vmapI)-Tin);
% Compute q and form differences at faces
q = phys.Dx.*( mesh.rx.*(line.Dr*var) ...
    - line.LIFT*(mesh.fScale.*(mesh.nx .* dT/2.0)) );
dq = (q(mesh.vmapM)-q(mesh.vmapP))/2.0;

dq(mapI) = 0.0;
% form field differences at faces
alpha=0;
maxvel = max(max(abs(a)));

dT2 = a(mesh.vmapM).*( var(mesh.vmapM)-var(mesh.vmapP) ); %.* ...
flux = mesh.nx.*(dT2./2.0 - dq) - maxvel/2.0.*dT;
%     (u(mesh.vmapM) .* mesh.nx-(1-alpha)*abs(u(mesh.vmapM) .* mesh.nx))/2;

dfdr = line.Dr*(phys.u.*var - q);

% compute right hand sides of the semi-discrete PDE
rhsu = -mesh.rx.*(dfdr) + line.LIFT*( mesh.fScale.*flux );

end% func
