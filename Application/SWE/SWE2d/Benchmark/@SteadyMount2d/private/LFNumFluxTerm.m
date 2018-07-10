function [ FluxS ] = LFNumFluxTerm( hmin, gra, fm, fp, nx, ny )
% normal flux
qnM = + fm(:, :, 2) .* nx + fm(:, :, 3) .* ny;
qvM = - fm(:, :, 2) .* ny + fm(:, :, 3) .* nx;
qnP = + fp(:, :, 2) .* nx + fp(:, :, 3) .* ny;
qvP = - fp(:, :, 2) .* ny + fp(:, :, 3) .* nx;

% evaluate local eigenvalues of the flux Jacobian
lamda = max( ...
    sqrt( gra * fm(:, :, 1) ) + sqrt( qnM.^2 + qvM.^2 )./fm(:, :, 1), ...
    sqrt( gra * fp(:, :, 1) ) + sqrt( qnP.^2 + qvP.^2 )./fp(:, :, 1) );

EM = NodalFluxTerm( gra, fm(:, :, 1), qnM, qvM );
EP = NodalFluxTerm( gra, fp(:, :, 1), qnP, qvP );

FluxS(:, :, 1) = ( EM(:, :, 1) + EP(:, :, 1) ) .* 0.5 ...
    - lamda .* ( fp(:, :, 1) - fm(:, :, 1) );
Fqn = ( EM(:, :, 2) + EP(:, :, 2) ) .* 0.5 - lamda .* ( qnP - qnM );
Fqv = ( EM(:, :, 3) + EP(:, :, 3) ) .* 0.5 - lamda .* ( qvP - qvM );

FluxS(:, :, 2) = (Fqn .* nx - Fqv .* ny);
FluxS(:, :, 3) = (Fqn .* ny + Fqv .* nx);

end

function E = NodalFluxTerm( gra, h, hu, hv )
E(:, :, 1) = hu;
E(:, :, 2) = 0.5 * gra * h.^2 + hu.^2 ./ h;
E(:, :, 3) = hu .* hv ./ h;
end
