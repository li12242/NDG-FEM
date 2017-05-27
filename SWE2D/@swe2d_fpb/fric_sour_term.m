function [ sf ] = fric_sour_term( obj, f_Q )
%FRIC_SOUR_TERM Summary of this function goes here
%   Detailed explanation goes here

sf(:,:,1) = zeros(obj.mesh.cell.Np, obj.mesh.K);
u = sqrt( f_Q(:, :, 2).^2 + f_Q(:, :, 3).^2 );
sf(:,:,2) = -obj.k.*f_Q(:,:,2).*u./f_Q(:,:,1).^2;
sf(:,:,3) = -obj.k.*f_Q(:,:,3).*u./f_Q(:,:,1).^2;

sf(:, ~obj.wetflag, 2) = 0;
sf(:, ~obj.wetflag, 3) = 0;
end

