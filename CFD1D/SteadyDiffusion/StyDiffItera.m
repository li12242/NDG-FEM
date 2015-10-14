function [c, q] = StyDiffItera(mesh, c, q)
% 1D Steady diffusion iteration

line = mesh.Shape;
vmapI = mesh.vmapM(mesh.mapI);
vmapO = mesh.vmapM(mesh.mapO);

qs = mesh.nx.*(q(mesh.vmapM) - ( c(mesh.vmapM) - c(mesh.vmapP) ));
qin = q(vmapI); cin = -c(vmapI);
qs(mesh.mapI) = mesh.nx(mesh.mapI).*(q(vmapI) - ( c(vmapI) - cin ));
qout =q(vmapO); cout = -c(vmapO);
qs(mesh.mapO) = mesh.nx(mesh.mapO).*(q(vmapO) - ( c(vmapO) - cout));
dq = mesh.nx.*(q(mesh.vmapM)) - qs;

q = - (line.Dr - diag(diag(line.Dr)))*q + (line.invM*line.FaceMassMatrixSmall*(mesh.fScale.*dq) - 1)./mesh.rx;
q = diag(diag(1./line.Dr))*q;
% boundary condition

c = -(line.Dr - diag(diag(line.Dr)))*c + q./mesh.rx;
c = diag(diag(1./line.Dr))*c;

end% func