function [Fh, Fq] = SWEFlux(h, q, hFlux)
% calculate flux terms of SWE

g = 9.81; isWet = (h > hFlux);
u = zeros(size(h)); 

u(isWet) = q(isWet)./h(isWet);
Fh = q;
Fq = g*h.^2./2 + u.^2.*h;

% for dry elements, no flow flux
Fq(~isWet) = 0.0;
end% func