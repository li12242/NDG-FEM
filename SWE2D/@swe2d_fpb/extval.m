function [ f_extQ ] = extval( obj, time )
%EXTVAL Summary of this function goes here
%   Detailed explanation goes here

gra = obj.gra; 
k = obj.k; % 线性摩阻系数
h0 = obj.h0; 
a = obj.a; 
B = obj.B;

p = sqrt(8*gra*h0./a^2);
s = sqrt(p^2 - k^2);

eta = h0 - B^2./2./gra*exp(-k*time) ...
    - B./gra*exp(-k*time/2)*(k/2*sin(s*time) + s*cos(s*time)).*obj.mesh.x ...
    - B./gra*exp(-k*time/2)*(k/2*cos(s*time) + s*sin(s*time)).*obj.mesh.y;

h = eta - obj.bot; h(h<0) = 0; f_extQ(:,:,1) = h;

u = B*exp(-k*time/2)*sin(s*time); f_extQ(:,:,2) = u.*h;
v = B*exp(-k*time/2)*sin(s*time); f_extQ(:,:,3) = v.*h;

end

