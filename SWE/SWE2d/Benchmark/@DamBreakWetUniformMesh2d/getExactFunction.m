function [ fext ] = getExactFunction( obj, time )
theta = obj.theta;
fext = cell( obj.Nmesh, 1 );

for m = 1:obj.Nmesh
    fext{m} = zeros(obj.meshUnion(m).cell.Np, obj.meshUnion(m).K, obj.Nfield);
    xc = obj.damPosition;
    yc = 0;
    dx = obj.meshUnion(m).x - xc;
    dy = obj.meshUnion(m).y - yc;
    x1d = xc + ( cos(-theta) * dx - sin(-theta) * dy );
    [ h, u ] = getExactDamBreak1d(obj, obj.meshUnion(m), x1d, time);
    hu = h.*u;
    fext{m}(:,:,1) = h;
    fext{m}(:,:,2) = cos(obj.theta) * hu;
    fext{m}(:,:,3) = sin(obj.theta) * hu;
end

end

%======================================================================
%> @brief Calculate the exact solution of wet dam break problem for 1d case
%>
%> More detailed description.
%>
%> @param arg1 First argument
%> @param arg2 Second argument
%>
%> @retval out1 return value for the first output variable
%> @retval out2 return value for the second output variable
%======================================================================
%> This function is part of the NDGOM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
%======================================================================
function [ h, u ] = getExactDamBreak1d(obj, mesh, x, time)
hr  = obj.h1;
hl  = obj.h0;
xc    = obj.damPosition;
theta = (x-xc)/time;
g     = obj.gra;
sghl  = sqrt(g*hl);

s     = ShockSpeed1d( g,hl,hr );
hm    = WaterHeightMiddle( g,hr,s );
um    = VelocityMiddle( g,hr,s );
sghm  = um - sqrt(g*hm);

h     = ones(mesh.cell.Np, mesh.K)*hr;
u     = zeros(mesh.cell.Np, mesh.K);
% wet part
ind    = theta < -sghl;
h(ind) = hl;
u(ind) = 0;
% rare wave
ind    = (-sghl<=theta) & (theta<sghm);
h(ind) = 1/9/g.*(2*sghl - theta(ind)).^2;
u(ind) = 2/3*(sghl + theta(ind));
% middle part
ind    = (theta >= sghm) & (theta < s);
h(ind) = hm;
u(ind) = 2/3*(sghl + sghm);
% right part
ind    = (theta >= s);
h(ind) = hr;
u(ind) = 0;
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
