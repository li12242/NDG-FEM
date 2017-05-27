function [ sf ] = fric_sour_term( obj, f_Q )
%FRIC_SOUR_TERM Summary of this function goes here
%   Detailed explanation goes here

sf(:,:,1) = zeros(obj.mesh.cell.Np, obj.mesh.K);
u = sqrt(f_Q(:,:,2).^2 + f_Q(:,:,3).^2);

sf(:,:,2) = -obj.tau.*u;
sf(:,:,3) = -obj.tau.*u;

sf(:, ~obj.wetflag, 2) = 0;
sf(:, ~obj.wetflag, 3) = 0;
end

