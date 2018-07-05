function frhs2d = EvaluatePCE2d_BoundaryKernel( obj, edge, fphys2d, fext2d )
%EVALUATEPCE2D_BOUNDARYKERNEL Summary of this function goes here
%   Detailed explanation goes here

[ fm, fp ] = edge.matEvaluateSurfValue( fphys2d );

% apply clamped boundary condition
ind = ( edge.ftype == enumBoundaryCondition.Clamped );
% fp(:, ind, 4) = fext2d(:, ind, 4); 
fp(:, ind, 1) = fext2d(:, ind, 1);

% apply slip wall boundary condition
ind = ( edge.ftype == enumBoundaryCondition.SlipWall );
Hun =  fm( :, ind, 2 ) .* edge.nx(:, ind) + fm( :, ind, 3).* edge.ny(:, ind);
Hvn = -fm( :, ind, 2 ) .* edge.ny(:, ind) + fm( :, ind, 3).* edge.nx(:, ind);

fp(:, ind, 2) = - Hun .* edge.nx(:, ind) - Hvn .* edge.ny(:, ind);
fp(:, ind, 3) = - Hun .* edge.ny(:, ind) + Hvn .* edge.nx(:, ind);

lambda = max( sqrt( obj.gra .* fm(:, :, 4) ), sqrt( obj.gra .* fp(:, :, 4) ) );

FluxM = fm(:, :, 2) .* edge.nx + fm(:, :, 3) .* edge.ny;
FluxP = fp(:, :, 2) .* edge.nx + fp(:, :, 3) .* edge.ny;
FluxS = 0.5 * ( FluxM + FluxP - lambda .* ( fp(:, :, 1) - fm(:, :, 1) ) );

frhs2d = edge.matEvaluateStrongFromEdgeRHS( FluxM, FluxS );
end

