function s = ShockSpeed1d( gra,hl,hr )
%SHOCKSPEED1D Determine the shock speed for dam break problem
%   Refer to Zhao, Guo (2014) for more details.

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

function hm = WaterHeightMiddle(gra, hr, s)
hm = hr./2.*( sqrt(1+8.*s.^2/gra/hr)-1 );
end% func

function um = VelocityMiddle(gra, hr, s)
um = s-gra.*hr./4./s.*( sqrt(1+8.*s.^2/gra/hr)+1 );
end% func

