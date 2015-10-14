function dc = LinConBC(mesh, c, dc, u, time)
% 1D Linear Convection Test Case
% Purpose  : Impliment Boundary Condition

% Inflow
cin = 0;

vmapI = mesh.vmapM(mesh.mapI);
cs = mesh.nx(mesh.mapI).*u.*(c( vmapI )+cin)./2 ...
    - abs(u*mesh.nx(mesh.mapI)).*(cin - c(vmapI))./2;

dc(mesh.mapI) = (mesh.nx(mesh.mapI).*u.*c(mesh.vmapM(mesh.mapI)) - u.*cs);
% dc(mesh.mapI) = 0.0;
% Outflow
% dc(mesh.mapO)
end% func