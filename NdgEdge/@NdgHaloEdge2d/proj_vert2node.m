function [ fnode ] = proj_vert2node( obj, fvert )
%PROJ_VERT2NODE Summary of this function goes here
%   Detailed explanation goes here

fnode = 0.5 * ( (1 - obj.r) * fvert(1,:) + (1 + obj.r) * fvert(2,:) );

end
