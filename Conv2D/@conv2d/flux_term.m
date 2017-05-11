function [ E, G ] = flux_term( obj, f_Q )
%FLUX_TERM Summary of this function goes here
%   Detailed explanation goes here

E = obj.u.*f_Q;
G = obj.v.*f_Q;
end

