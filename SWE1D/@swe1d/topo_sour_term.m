function [ sb ] = topo_sour_term( obj, f_Q )
%SOURCE_TERM µ×ÆÂÔ´Ïî
%   Detailed explanation goes here

sb(:,:,1) = zeros(obj.mesh.cell.Np, obj.mesh.K);
sb(:,:,2) = -obj.gra.*(f_Q(:,:,1) + obj.bot)...
    .*obj.mesh.rx.*(obj.mesh.cell.Dr*(obj.bot));

sb(:, ~obj.wetflag, 2) = 0;
end

