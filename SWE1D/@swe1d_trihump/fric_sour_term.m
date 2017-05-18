function [ sf ] = fric_sour_term( obj, f_Q )
%FRIC_SOUR_TERM Summary of this function goes here
%   Detailed explanation goes here

sf(:,:,1) = zeros(obj.mesh.cell.Np, obj.mesh.K);
sf(:,:,2) = -obj.gra.*obj.n.*f_Q(:,:,2).*abs(f_Q(:,:,2))./f_Q(:,:,1).^(7/3);

sf(:, ~obj.wetflag, 2) = 0;
end

