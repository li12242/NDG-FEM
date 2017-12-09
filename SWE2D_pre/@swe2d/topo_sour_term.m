function [ sb ] = topo_sour_term( obj, f_Q )
%TOPO_SOUR_TERM Summary of this function goes here
%   Detailed explanation goes here

sb = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);

sb(:, :, 2) = -obj.gra.*( f_Q(:,:,1)+obj.bot ).*obj.bx;
sb(:, :, 3) = -obj.gra.*( f_Q(:,:,1)+obj.bot ).*obj.by;

sb(:, ~obj.wetflag, 2) = 0;
sb(:, ~obj.wetflag, 3) = 0;
end
