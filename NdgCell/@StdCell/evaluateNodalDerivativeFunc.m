function [ fDr, fDs, fDt ] = evaluateNodalDerivativeFunc( obj, r, s, t )
Nr = numel( r );
Vr = zeros(Nr, obj.Np);
Vs = zeros(Nr, obj.Np);
Vt = zeros(Nr, obj.Np);
for n = 1:obj.Np
    [Vr(:, n), Vs(:, n), Vt(:, n)] = ...
        obj.evaluateDerivativeOrthogonalFunc(obj.N, n, r, s, t);
end
fDr = Vr/obj.V;
fDs = Vs/obj.V;
fDt = Vt/obj.V;
end