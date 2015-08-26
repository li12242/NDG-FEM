function [rhsQ, lambda] = BurgersRHS(mesh, Q)
% 1D Burgers¡¯ equation
% right hand side

line = mesh.Shape;

QM = Q(mesh.vmapM); QP = Q(mesh.vmapP);
lambda = 0.5.*(QM+QP).*mesh.nx;

FM = BurgersFlux(QM); FP = BurgersFlux(QP);
Fs = 0.5*(1+sign(lambda)).*FM + 0.5*(1-sign(lambda)).*FP;
dF = mesh.nx.*(FM - Fs);

F = BurgersFlux(Q);

rhsQ= -mesh.rx.*(line.Dr*F)+line.invM*line.FaceMassMatrixSmall*(mesh.fScale.*dF);
end

function F = BurgersFlux(Q)
F = 0.5*Q.^2;
end