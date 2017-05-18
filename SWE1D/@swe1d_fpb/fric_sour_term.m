function [ sf ] = fric_sour_term( obj, f_Q )
%FRIC_SOUR_TERM Summary of this function goes here
%   Detailed explanation goes here

sf(:,:,1) = zeros(obj.mesh.cell.Np, obj.mesh.K);
sf(:,:,2) = -obj.tau.*f_Q(:,:,2);

sf(:, ~obj.wetflag, 2) = 0;
end

