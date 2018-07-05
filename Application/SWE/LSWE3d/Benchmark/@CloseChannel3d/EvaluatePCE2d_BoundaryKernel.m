function frhs2d = EvaluatePCE2d_BoundaryKernel( obj, edge, fphys2d, fext2d )
%EVALUATEPCE2D_BOUNDARYKERNEL Summary of this function goes here
%   Detailed explanation goes here

[ fm, fp ] = edge.matEvaluateSurfValue( fphys2d );

% apply slip wall boundary condition
Hun =  fm( :, :, 2 ) .* edge.nx + fm( :, :, 3).* edge.ny;
Hvn = -fm( :, :, 2 ) .* edge.ny + fm( :, :, 3).* edge.nx;

fp(:, :, 2) = - Hun .* edge.nx - Hvn .* edge.ny;
fp(:, :, 3) = - Hun .* edge.ny + Hvn .* edge.nx;

FluxM = fm(:, :, 2) .* edge.nx + fm(:, :, 3) .* edge.ny;
FluxP = fp(:, :, 2) .* edge.nx + fp(:, :, 3) .* edge.ny;
FluxS = 0.5 * ( FluxM + FluxP );

frhs2d = edge.matEvaluateStrongFromEdgeRHS( FluxM, FluxS );
end

