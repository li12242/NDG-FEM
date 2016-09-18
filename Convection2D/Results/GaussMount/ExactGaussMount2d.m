function [ var ] = ExactGaussMount2d( x, y, t )
%EXACTGAUSSMOUNT2D Exact solution of GaussMount test case
%   Detailed explanation goes here
var = zeros(size(x));
r0  = 0.15;
x0  = 0.25; 
y0  = 0.5;
r2  = sqrt((x-x0).^2+(y-y0).^2)./r0;
ind = ( r2<=1.0);
var(ind) = (1+cos(r2(ind)*pi))./2;
end