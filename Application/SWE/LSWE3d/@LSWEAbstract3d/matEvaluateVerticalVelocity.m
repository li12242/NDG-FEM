function [ W ] = matEvaluateVerticalVelocity( obj, mesh3d, fphys2d, fphys3d )
%MATEVALUATEVERTICALVELOCITY Summary of this function goes here
%   Detailed explanation goes here

% fphys3d(:, :, 6) = mesh3d.Extend2dField( fphys2d(:, :, 1) - fphys2d(:, :, 5) );

% HUt = mesh3d.Extend2dField( fphys2d(:, :, 2) ) - fphys3d(:, :, 5) .* fphys3d(:, :, 1);
% HVt = mesh3d.Extend2dField( fphys2d(:, :, 3) ) - fphys3d(:, :, 5) .* fphys3d(:, :, 2);

HUt = - fphys3d(:, :, 6) .* fphys3d(:, :, 1);
HVt = - fphys3d(:, :, 6) .* fphys3d(:, :, 2);

dHUdx = mesh3d.rx .* ( mesh3d.cell.Dr * HUt ) + mesh3d.sx .* ( mesh3d.cell.Ds * HUt );
dHVdy = mesh3d.ry .* ( mesh3d.cell.Dr * HVt ) + mesh3d.sy .* ( mesh3d.cell.Ds * HVt );

W = mesh3d.VerticalIntegralField( dHUdx + dHVdy );
end

