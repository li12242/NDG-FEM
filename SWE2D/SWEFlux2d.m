function [Fx, Fy] = SWE2DFlux(Q)
% two-dimension shallow water flux
g = 9.8; Fx = zeros(size(Q)); Fy = zeros(size(Q));

h = Q(:,:,1); hu = Q(:,:,2); hv = Q(:,:,3);

% case h=0
iswet = (h(:)>10^-10); temp = zeros(size(h));
Fx(:,:,1) = hu;
temp(iswet) = hu(iswet).^2./h(iswet) + 0.5.*g.*h(iswet).^2; 
Fx(:,:,2) = temp; temp(iswet) = 0;
temp(iswet) = hu(iswet).*hv(iswet)./h(iswet);
Fx(:,:,3) = temp; temp(iswet) = 0;

Fy(:,:,1) = hv; 
temp(iswet) = hu(iswet).*hv(iswet)./h(iswet); 
Fy(:,:,2) = temp; temp(iswet) = 0;
temp(iswet) = hv(iswet).^2./h(iswet) + 0.5.*g.*h(iswet).^2;
Fy(:,:,3) = temp;
end