function [ var ] = ExactSquareMount2d( x,y,t )
%EXACTSQUAREMOUNT2D Exact solution of square mount test case
%   Detailed explanation goes here

var = zeros(size(x));
b   = 1/12;
x0  = 0.50; xc = mean(x); 
y0  = 0.75; yc = mean(y);

flag = ( abs(xc - x0)<b & abs(yc - y0)<b );
var(:, flag) = 1;
end

