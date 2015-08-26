function F = SWEFlux(Q)
% F = SWEFlux(Q);
% calculate the flux of shallow water equation
g = 9.8;
h = Q(:,:,1); hu = Q(:,:,2); u = hu./h;
F(:,:,1) = hu; F(:,:,2) = hu.*u + 1/2*g*h.^2;
end