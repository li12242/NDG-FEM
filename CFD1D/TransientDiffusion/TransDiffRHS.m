function rhs = TransDiffRHS(mesh, c, time)
% 1D Transient Diffusion Test Case
% Purpose: Evaluate RHS flux in strong form

line = mesh.Shape;

% source 
% Q = ones(size(c));
vmapI = mesh.vmapM(mesh.mapI);
vmapO = mesh.vmapM(mesh.mapO);

cs = mesh.nx.*(c(mesh.vmapM) + c(mesh.vmapP))./2;
cin = -c(vmapI);    cs(mesh.mapI) = mesh.nx(mesh.mapI)*(cin + c(vmapI))./2;
cout = -c(vmapO);   cs(mesh.mapO) = mesh.nx(mesh.mapO)*(cout + c(vmapO))./2;
dc = mesh.nx.*(c(mesh.vmapM)) - cs;
q = -mesh.rx.*(line.Dr*c) + line.invM*line.FaceMassMatrixSmall*(mesh.fScale.*dc);

qs = mesh.nx.*(q(mesh.vmapM) + q(mesh.vmapP))./2;
qin = q(vmapI);     qs(mesh.mapI) = mesh.nx(mesh.mapI)*(qin + q(vmapI))./2;
qout =q(vmapO);     qs(mesh.mapO) = mesh.nx(mesh.mapO)*(qout +q(vmapO))./2;
dq = mesh.nx.*(q(mesh.vmapM)) - qs;
rhs = -mesh.rx.*(line.Dr*q) ...
    + line.invM*line.FaceMassMatrixSmall*(mesh.fScale.*dq) + 1;
end% func