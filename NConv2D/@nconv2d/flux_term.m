function [ E, G ] = flux_term( obj, f_Q )
%FLUX_TERM Summary of this function goes here
%   Detailed explanation goes here

E = f_Q.^2 * 0.5;
G = E;
end

