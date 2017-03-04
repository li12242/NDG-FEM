function [h, q, u] = exactParabolicBowlSolution(t, x, bedElva)
% analysitical solution of parabolic bowl solution

g = 9.8; B = 5; h0 = 10; a = 3000;
w = sqrt(2*g*h0)./a;

z = (-B^2*cos(2*w*t) - B^2 -(4*B*w)*cos(w*t).*x)./(4*g);
h = z - bedElva;
h(h<0) = 0;

u = zeros(size(x));
u(h>0) = B*a*w./sqrt(2*h0*g)*sin(w*t);
u(h<=0) = 0;

q = u.*h;
end% func