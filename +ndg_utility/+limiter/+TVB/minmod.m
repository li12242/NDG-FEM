function [ m ] = minmod( v )
%MINMOD Summary of this function goes here
%   Detailed explanation goes here

num = size(v,1);
m = zeros(1, size(v,2));
s = sum(sign(v), 1)/num ;

ids = ( abs(abs(s)-1)<1e-10 );
tmp = s.*min( abs(v), [], 1);
m(ids) = tmp(ids);
end

