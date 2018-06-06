function E = flux_term( obj, f_Q )
%FLUX_TERM Summary of this function goes here
%   Detailed explanation goes here

E = f_Q .* obj.u;
end