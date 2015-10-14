function [Fh, Fq] = SWEFlux(h, q)

g = 9.8; hDelta = 10^-6; isWet = h>hDelta;
u = zeros(size(h)); 

u(isWet) = q(isWet)./h(isWet);
Fh = q;
Fq = g*h.^2./2 + u.^2.*h;
end% func