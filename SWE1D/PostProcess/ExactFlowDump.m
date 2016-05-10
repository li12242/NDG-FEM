function [x, z, h] = ExactFlowDump(condition)
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
        q0 = 0.18; h0 = 0.413203097187649; 
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

% 
% flag = z > eps;
% C = -q0^2/(2*g*h0^2) - h0;
% b = z(1);
% str = ['h^3 + ', num2str(b+C), '*h^2 + ', num2str(q0^2/2/g)];
% h1 = solve(str, 'h'); h1 = double(h1);
% [~, I] = min(abs((h1 - h0)));
% h(~flag) = h1(I);
% 
% % none zero bottom
% nodeIndex = find(flag);
% for i = 1:numel(nodeIndex)
%     b = z(nodeIndex(i));
%     str = ['h^3 + ', num2str(b+C), '*h^2 + ', num2str(q0^2/2/g)];
%     h1 = solve(str, 'h'); h1 = double(h1);
%     [~, I] = min(abs((h1 - h0)));
%     I = 1;
%     h(nodeIndex) =  h1(I);
% end

plot(x, z, 'k', x, z+h, 'b-')