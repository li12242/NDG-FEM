function [Fh, Fq] = SWEFlux(phys, h, q, isWet)
% calculate flux terms of SWE

%% Parameters
minDepth = phys.minDepth;
g        = phys.gra;
%% Wet/Dry status
wetNode  = (h > minDepth);
%% flow rate
u        = zeros(size(h)); 
u(wetNode) = q(wetNode)./h(wetNode);

%% Flux terms
Fh = q;
Fq = g*h.^2./2 + u.^2.*h;
% for dry elements, no flow flux
Fq(:, ~isWet) = 0.0;
end% func