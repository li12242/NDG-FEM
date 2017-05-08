function [ r,s,t ] = node_coor_func( obj, order )
%GETCOOR Summary of this function goes here
%   Detailed explanation goes here

np = order+1;
[x,~] = Polylib.zwglj(np);

% 首先按照x坐标循环，然后按照y坐标循环
r = x*ones(1, np);
s = ones(np, 1)*x';

r = r(:); s = s(:);
t = zeros(size(r));
end

