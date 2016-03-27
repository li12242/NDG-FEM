function [x, h] = ExactFlowDump(condition)
% exact solution of flow dump
% Input:
%   condition - [1]: Subcritical flow
%               [2]: Supercritical flow
%               [3]: Transcritical flow

x1 = 0; x2 = 25; % Flow over dump

ne = 500; np = ne + 1;
x = linspace(x1, x2, np);

% bottom topography
flag = (x >= 8) & (x <=12);
z = zeros(size(x));
z(flag) = 0.2 - 0.05*(x(flag) -10).^2;

% inflow condition
switch condition
    case 1
        q0 = 0.18; h0 = 0.5; 
    case 2
        q0 = 25.0567; h0 = 2.0; 
    case 3
        q0 = 0.18; h0 = 0.33; 
end
g = 9.8;

% 
h = zeros(size(x));
h(1) = h0;
for in = 2:np
    coef = ( 1 - q0^2./(g*h(in-1).^3) );
    xdelta = x(in) - x(in - 1);
    h(in) = h(in - 1) + xdelta./coef*( -(z(in) - z(in-1))./xdelta );
end% for

plot(x, z, 'k', x, z+h, 'b-')