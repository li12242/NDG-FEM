function rhs = rhs_term(obj, c, time)
%RHS_TERM Summary of this function goes here
%   Detailed explanation goes here

% volume integral
E = flux_term(obj, obj.u, c);
rhs = -obj.mesh.rx.*(obj.mesh.cell.Dr*E);

% surface integral
eidM = obj.mesh.eidM;
eidP = obj.mesh.eidP;
cM = c(eidM);
cP = c(eidP);
uM = obj.u(eidM);
uP = obj.u(eidP);

cP = nei_node_val(obj, cM, cP, obj.c_ext(obj.mesh.eidM), obj.mesh.eidtype);
fluxS = num_flux(obj, cM, uM, cP, uP, obj.mesh.nx); % numerical flux

fluxM = flux_term(obj, uM, cM); % local flux term
dflux = fluxM.*obj.mesh.nx - fluxS; % deviation of flux term

rhs = rhs + obj.mesh.cell.LIFT*( obj.mesh.eidfscal.*(dflux) );
end

