function [ h ] = DamBreakDry1d_ExactH( x,t )
%DAMBREAK1D_EXACTH exact solution of one dimensional DamBreak problem
%   Detailed explanation goes here

h0    = 10;
xc    = 500;
theta = (x - xc)/t;
g     = 9.81; 
sgh   = sqrt(g*h0);

h     = zeros(size(x));

% wet part
ind    = theta < -sgh;
h(ind) = h0;

ind    = theta > 2*sgh;
h(ind) = 0;

ind    = (theta >= -sgh) & (theta <= 2*sgh);
h(ind) = 1/9/g*(theta(ind) - 2*sgh).^2;
end

