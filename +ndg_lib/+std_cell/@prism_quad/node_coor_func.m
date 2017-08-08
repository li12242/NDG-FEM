function [ Np, r, s, t ] = node_coor_func( obj, N )
%NODE_COOR_FUNC Summary of this function goes here
%   Detailed explanation goes here

Np = (N+1)^3;

[ Nquad, qr, qs, qt ] = quad_node_coor( N );
[ Nline, lr, ls, lt ] = line_node_coor( N );

r = qr * ones(1, Nline);
s = qs * ones(1, Nline);
t = ones( Nquad, 1 ) * lr';

r = r(:);
s = s(:);
t = t(:);
end

function [ Np,r,s,t ] = quad_node_coor(N)
Np = N+1;
[x,~] = Polylib.zwglj(Np);

% 首先按照x坐标循环，然后按照y坐标循环
r = x*ones(1, Np);
s = ones(Np, 1)*x';

r = r(:); s = s(:); 
t = zeros(size(r));
Np = Np * Np;
end

function [ Np, r, s, t ] = line_node_coor( N )
Np = N+1;
[r,~] = Polylib.zwglj(Np);
s = zeros(Np, 1);
t = zeros(Np, 1);
end
