function [ sb ] = topo_sour_term( obj, f_Q )
%TOPO_SOUR_TERM Summary of this function goes here
%   Detailed explanation goes here

sb = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);



sb(:, :, 2) = -obj.gra.*( f_Q(:,:,1)+obj.bot ).*( ...
    obj.mesh.rx.*(obj.mesh.cell.Dr*obj.bot) + ...
    obj.mesh.sx.*(obj.mesh.cell.Ds*obj.bot) );
sb(:, :, 3) = -obj.gra.*( f_Q(:,:,1)+obj.bot ).*( ...
    obj.mesh.ry.*(obj.mesh.cell.Dr*obj.bot) + ...
    obj.mesh.sy.*(obj.mesh.cell.Ds*obj.bot) );

end

