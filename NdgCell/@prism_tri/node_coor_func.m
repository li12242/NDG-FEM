function [ Np,r,s,t ] = node_coor_func( obj, N )
%NODE_COOR_FUNC Summary of this function goes here
%   Detailed explanation goes here

Np = (N+1) * (N+1)*(N+2)/2;
% calculate the triangle
[ Ntri, tr, ts, tt ] = tri_node_coor( N );
[ Nlin, lr, ls, lt ] = line_node_coor( N );

r = tr * ones(1, Nlin);
s = ts * ones(1, Nlin);
t = ones( Ntri, 1 ) * lr';

r = r(:);
s = s(:);
t = t(:);
end

function [ Np, r, s, t ] = line_node_coor( N )
Np = N+1;
[r,~] = Polylib.zwglj(Np);
s = zeros(Np, 1);
t = zeros(Np, 1);
end



