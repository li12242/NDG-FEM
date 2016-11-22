function [x, y, rx, sx, ry, sy, J] = GetEleGeometric( obj, VX, VY )
%GETELEGEOMETRIC Summary of this function goes here
%   Detailed explanation goes here

% get element Geometric Factor
% Input:    vx - Vertic Coordinate, size [3(nVertice) x nElement]
%           vy - Vertic Coordinate, size [3(nVertice) x nElement]
% Output:   x  - node coordinate
%           rx - dr/dx at nodes
%           J  - jacobi factor
assert((size(VX,1)==4), 'Quad GetEleGeometric: input vertex VX has fault')
assert((size(VY,1)==4), 'Quad GetEleGeometric: input vertex VY has fault')

% coordinate mapping
x = (1-obj.r).*(1-obj.s)./4*VX(1,:) + (1+obj.r).*(1-obj.s)./4*VX(2,:)...
    +(1+obj.r).*(1+obj.s)./4*VX(3,:) +(1-obj.r).*(1+obj.s)./4*VX(4,:);
y = (1-obj.r).*(1-obj.s)./4*VY(1,:) + (1+obj.r).*(1-obj.s)./4*VY(2,:)...
    +(1+obj.r).*(1+obj.s)./4*VY(3,:) +(1-obj.r).*(1+obj.s)./4*VY(4,:);

xr = obj.Dr*x; xs = obj.Ds*x; yr = obj.Dr*y; ys = obj.Ds*y; 
J = -xs.*yr + xr.*ys;
rx = ys./J; sx =-yr./J; ry =-xs./J; sy = xr./J;

end

