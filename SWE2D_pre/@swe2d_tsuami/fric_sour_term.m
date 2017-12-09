function [ sf ] = fric_sour_term( obj, f_Q )
%FRIC_SOUR_TERM Calculate the friction source term.
%   The friction source term is $S = [0, -ghS_{fx}, -ghS_{fy}]$, where
%   $S_{fx}$ and $S_{fy}$ are the friction losses, and are estimated by 
%   the Manning law given by
%
%   $$S = [0, -ghS_{fx}, -ghS_{fy}]^T$$
%
%    (Eskilsson & Sherwin, 2004)

sf(:,:,1) = zeros(obj.mesh.cell.Np, obj.mesh.K);
qn  = sqrt( f_Q(:,:,2).^2+f_Q(:,:,3).^2 );
sf(:,:,2) = -obj.gra.*obj.n.*f_Q(:,:,2).*qn./( f_Q(:,:,1).^(7/3) );
sf(:,:,3) = -obj.gra.*obj.n.*f_Q(:,:,3).*qn./( f_Q(:,:,1).^(7/3) );

sf(:, ~obj.wetflag, 2) = 0;
sf(:, ~obj.wetflag, 3) = 0;
end

