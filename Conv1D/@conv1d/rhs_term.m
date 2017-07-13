function rhs = rhs_term( obj, f_Q )
%RHS_TERM Summary of this function goes here
%   Detailed explanation goes here

E = flux_term(obj, f_Q); % volume flux term
dflux = surf_term( obj, f_Q ); % surface flux deviation

rhs = -obj.mesh.rx.*(obj.mesh.cell.Dr*E) + ...
    obj.mesh.cell.LIFT*( obj.mesh.Js.*dflux )./obj.mesh.J;

end

