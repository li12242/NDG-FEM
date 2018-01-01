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
    [ h, u ] = getExtFunc1d(obj, obj.meshUnion(m), x1d, time);
    hu = h.*u;
    fext{m}(:,:,1) = h;
    fext{m}(:,:,2) = cos(obj.theta) * hu;
    fext{m}(:,:,3) = sin(obj.theta) * hu;
end

end

function [h, u] = getExtFunc1d(obj, mesh, x, time)
h = zeros(size(x));
u = zeros(size(x));
temp = (x - obj.damPosition)/time;
c0 = sqrt(obj.gra .* obj.h0);
% left part
ind = (temp < -c0 );
h(ind) = obj.h0; 
u(ind) = 0;
% middle part
ind = (temp > -c0 ) & ( temp < 2*c0 );
h(ind) = (2*c0 - temp(ind) ).^2/obj.gra/obj.gra;
u(ind) = 2/3*(c0 + temp(ind) );

end