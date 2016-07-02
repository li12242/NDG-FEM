function [rhsH, rhsQ] = SWERHS(phys, mesh, h, q, bot, isWet)
% Righ Hand Side

line = mesh.Shape;
% numel flux
[Fhs, Fqs, ~] = SWEHLL(phys, mesh, h, q, isWet);
% source term
[Sh, Sq]      = SWESource(phys, mesh, bot, h, isWet);

% flux term
[Fh, Fq] = SWEFlux(phys, h, q, isWet);

% rhs
rhsH = mesh.rx.*(line.invM*(line.Dr'*(line.M*Fh))) ...
    - line.LIFT*(Fhs.*mesh.fScale) + Sh;
rhsQ = mesh.rx.*(line.invM*(line.Dr'*(line.M*Fq))) ...
    - line.LIFT*(Fqs.*mesh.fScale) + Sq;
end% func

function [Sh, Sq] = SWESource(phys, mesh, bot, h, isWet)

g    = phys.gra; 
line = mesh.Shape;
Sh   = zeros(size(h));
Sq   = -g.*h.*mesh.rx.*(line.Dr*bot);
Sq(:, ~isWet) = 0.0;
end% func



