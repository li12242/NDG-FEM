function [h, qx, qy] = ParabolicBowlExtSol(x, y, bot, t)
% Parameters
g     = 9.81;
alpha = 1.6*1e-7;
w     = sqrt(8*g*alpha);
X     = 1;
Y     = -0.41884;

% get exact solution
r2    = x.^2 + y.^2;
h     = 1/(X+Y*cos(w*t)) + alpha*(Y^2 - X^2).*r2/(X+Y*cos(w*t))^2 - bot;
h(h<0)=0;
u     = -(Y*w*sin(w*t))./(X+Y*cos(w*t)).*x/2;
v     = -(Y*w*sin(w*t))./(X+Y*cos(w*t)).*y/2;
qx    = u.*h;
qy    = v.*h;
end% func