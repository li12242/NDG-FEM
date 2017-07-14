function [ h ] = DamBreakWet2d_H( x,y,t )
%DAMBREAKWET2D_H Summary of this function goes here
%   Detailed explanation goes here

hr  = 2;
hl  = 10;
xc    = 500;
theta = (x-xc)/t;
g     = 9.81;
sghl  = sqrt(g*hl);

s     = ShockSpeed1d( g,hl,hr );
hm    = WaterHeightMiddle( g,hr,s );
um    = VelocityMiddle( g,hr,s );
sghm  = um - sqrt(g*hm);

h     = zeros(size(x));
% wet part
ind    = theta < -sghl;
h(ind) = hl;
% rare wave
ind    = (-sghl<=theta) & (theta<sghm);
h(ind) = 1/9/g.*(2*sghl - theta(ind)).^2;
% middle part
ind    = (theta>=sghm) & (theta<s);
h(ind) = hm;
% right part
ind    = (theta>=s);
h(ind) = hr;
end

function hm = WaterHeightMiddle(gra, hr, s)
hm = hr./2.*( sqrt(1+8.*s.^2/gra/hr)-1 );
end% func

function um = VelocityMiddle(gra, hr, s)
um = s-gra.*hr./4./s.*( sqrt(1+8.*s.^2/gra/hr)+1 );
end% func

