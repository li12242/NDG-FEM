function [ f_ext ] = ext_func( obj, t )
%EXT_FUNC The exact solutions for the dam break (dry) problem.
%   The exact solution is obtained from 
%   [1] Jakeman J. On numerical solutions of the shallow water wave
%       equations. 2006. P36.
%
A = [cos(-obj.theta), sin(-obj.theta); 
    -sin(-obj.theta), cos(-obj.theta)];
x = zeros(obj.mesh.cell.Np, obj.mesh.K);
temp = A*[ obj.mesh.x(:)'; obj.mesh.y(:)'; ];
x(:) = temp(1, :);

f_ext = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);
hr  = obj.h1;
hl  = obj.h0;
xc    = obj.dam_pos;
theta = (x-xc)/t;
g     = obj.gra;
sghl  = sqrt(g*hl);

s     = ShockSpeed1d( g,hl,hr );
hm    = WaterHeightMiddle( g,hr,s );
um    = VelocityMiddle( g,hr,s );
sghm  = um - sqrt(g*hm);

h     = zeros(obj.mesh.cell.Np, obj.mesh.K);
u     = zeros(obj.mesh.cell.Np, obj.mesh.K);
% wet part
ind    = theta < -sghl;
h(ind) = hl;
u(ind) = 0;
% rare wave
ind    = (-sghl<=theta) & (theta<sghm);
h(ind) = 1/9/g.*(2*sghl - theta(ind)).^2;
u(ind) = 2/3*(sghl + theta(ind));
% middle part
ind    = (theta>=sghm) & (theta<s);
h(ind) = hm;
u(ind) = 2/3*(sghl + sghm);
% right part
ind    = (theta>=s);
h(ind) = hr;
u(ind) = 0;
f_ext(:,:,1) = h; hu = u.*h;
f_ext(:,:,2) = cos(obj.theta) * hu;
f_ext(:,:,3) = -sin(obj.theta) * hu;
end

function hm = WaterHeightMiddle(gra, hr, s)
hm = hr./2.*( sqrt(1+8.*s.^2/gra/hr)-1 );
end% func

function um = VelocityMiddle(gra, hr, s)
um = s-gra.*hr./4./s.*( sqrt(1+8.*s.^2/gra/hr)+1 );
end% func

function s = ShockSpeed1d( gra,hl,hr )
%SHOCKSPEED1D Determine the shock speed for dam break (dry) problem
%   Refer to Zhao and Guo et al. (2014) for more details.
%
%   [1] Zhao L, Guo B, et al. International Journal for Numerical Methods
%       in Fluids 2014;75:815?34. 

TOLERR = 1e-12;
s1  = 1e-12; 
s2  = gra*(hr+hl);
BOUNDERR = abs(s2 - s1);

s   = (s1+s2)./2;
while(BOUNDERR>TOLERR)
    y  = Func(gra, hl, hr, s);
    y1 = Func(gra, hl, hr, s1);
    y2 = Func(gra, hl, hr, s2);
    % new bound
    if (y*y1>=0)
        s1 = s;
    elseif(y*y2>=0)
        s2 = s;
    else
        error('result %f=y(%f) is out of bound, [%f, %f]',y, s, y1, y2);
    end% if
    s = (s1+s2)./2;
    BOUNDERR = abs(s2 - s1);
    % 
end% while
end

function y = Func(gra, hl, hr, s)
um = VelocityMiddle(gra, hr, s);
hm = WaterHeightMiddle(gra, hr, s);
y  = um+2*sqrt(gra*hm)-2*sqrt(gra*hl);
end% func




