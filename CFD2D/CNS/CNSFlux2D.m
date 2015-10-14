function [F,G,rho,rhou,rhov,Ener] = CNSFlux2D(Q, gamma)

rho = Q(:,:,1); rhou = Q(:,:,2); rhov = Q(:,:,3); Ener = Q(:,:,4);
u = rhou./rho; v = rhov./rho; Pr = (gamma-1)*(Ener-0.5*rho.*(u.^2+v.^2));  

F = zeros(size(Q));
F(:,:,1) = rhou; F(:,:,2) = rhou.*u + Pr; F(:,:,3) = rhov.*u; F(:,:,4) = u.*(Ener + Pr);

G = zeros(size(Q));
G(:,:,1) = rhov; G(:,:,2) = rhou.*v; G(:,:,3) = rhov.*v + Pr; G(:,:,4) = v.*(Ener +Pr);
end