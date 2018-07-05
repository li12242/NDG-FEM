function frhs3d = Evaluate3d_HorizontalBoundaryKernel( obj, edge, fphys3d, fext3d )
%EVALUATE3D_HORIZONTALBOUNDARYKERNEL Summary of this function goes here
%   Detailed explanation goes here

[ fm, fp ] = edge.matEvaluateSurfValue( fphys3d );

un =  fm( :, :, 1 ) .* edge.nx + fm( :, :, 2).* edge.ny;
vn = -fm( :, :, 1 ) .* edge.ny + fm( :, :, 2).* edge.nx;

fp(:, :, 1) = - un .* edge.nx - vn .* edge.ny;
fp(:, :, 2) = - un .* edge.ny + vn .* edge.nx;

FluxM(:, :, 1) = obj.gra .* fm( :, :, 7 ) .* edge.nx;
FluxP(:, :, 1) = obj.gra .* fp( :, :, 7 ) .* edge.nx;
FluxM(:, :, 2) = obj.gra .* fm( :, :, 7 ) .* edge.ny;
FluxP(:, :, 2) = obj.gra .* fp( :, :, 7 ) .* edge.ny;

lambda = sqrt( obj.gra .* fm(:, :, 5) );

FluxS(:, :, 1) = 0.5 .* ( FluxM(:, :, 1) + FluxP(:, :, 1) - ...
    lambda .* ( fp( :, :, 1 ) - fm( :, :, 1 ) ) );
FluxS(:, :, 2) = 0.5 .* ( FluxM(:, :, 2) + FluxP(:, :, 2) - ...
    lambda .* ( fp( :, :, 2 ) - fm( :, :, 2 ) ) );

frhs3d = edge.matEvaluateStrongFormEdgeRHS( FluxM, FluxS );

end

