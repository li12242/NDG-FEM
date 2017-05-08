function [ var ] = ExactSolidBody2d( x, y, t )
%EXACTSOLIDBODY2D Exact solution of SolidBody test case
%   Detailed explanation goes here

var = zeros(size(x));
% slotted cylinder
x0  = 0.5; 
y0  = 0.75; 
r0  = 0.15;
r2  = sqrt((x-x0).^2+(y-y0).^2)./r0;
ind = ( r2<=1.0);
ind = ind & ((abs(x - x0)>=0.025) | (y >= 0.85));
var(ind) = 1.0;
% cone
x0  = 0.5; 
y0  = 0.25;
r2  = sqrt((x-x0).^2+(y-y0).^2)./r0;
ind = ( r2<=1.0);
var(ind) = 1-r2(ind);
% hump
x0  = 0.25;
y0  = 0.5;
r2  = sqrt((x-x0).^2+(y-y0).^2)./r0;
ind = ( r2<=1.0);
var(ind) = (1+cos(r2(ind)*pi))./4;
end

