function [ he ] = ObliqueHydraulicJump_H( x,y,t )
%OBLIQUEHYDRAULICJUMP_H Summary of this function goes here
%   Detailed explanation goes here

he = ones(size(x));

k = -tan(30/180*pi);

ind = y >= ( 30+k*(x-10) );
he(ind)  = 1.5049;
end

