function [Fh, Fq] = SWEFlux(h, q, isWet, physics)
% calculate flux terms of SWE


hPositive = physics.getVal('minDepth');
g = physics.getVal('gravity');

wetNode = (h > hPositive);
u = zeros(size(h)); 

u(wetNode) = q(wetNode)./h(wetNode);
Fh = q;
Fq = g*h.^2./2 + u.^2.*h;

% for dry elements, no flow flux
Fq(:, ~isWet) = 0.0;
end% func