function frhs3d = Evaluate3d_HorizontalBoundaryKernel( obj, edge, fphys3d, fext3d )
%EVALUATE3D_HORIZONTALBOUNDARYKERNEL Summary of this function goes here
%   Detailed explanation goes here

[ fm, fp ] = edge.matEvaluateSurfValue( fphys3d );

% apply clamped boundary condition
ind = ( edge.ftype == enumBoundaryCondition.Clamped );
fp(:, ind, 1) = fext3d(:, ind, 1);
fp(:, ind, 2) = fext3d(:, ind, 2);
fp(:, ind, 7) = fext3d(:, ind, 7);
% fp(:, ind, 6) = fext3d(:, ind, 6);

% apply slip wall boundary condition
ind = ( edge.ftype == enumBoundaryCondition.SlipWall );
Hun =  fm( :, ind, 1 ) .* edge.nx(:, ind) + fm( :, ind, 2).* edge.ny(:, ind);
Hvn = -fm( :, ind, 1 ) .* edge.ny(:, ind) + fm( :, ind, 2).* edge.nx(:, ind);

fp(:, ind, 1) = - Hun .* edge.nx(:, ind) - Hvn .* edge.ny(:, ind);
fp(:, ind, 2) = - Hun .* edge.ny(:, ind) + Hvn .* edge.nx(:, ind);

FluxM(:, :, 1) = obj.gra .* fm( :, :, 6 ) .* edge.nx;
FluxP(:, :, 1) = obj.gra .* fp( :, :, 6 ) .* edge.nx;
FluxM(:, :, 2) = obj.gra .* fm( :, :, 6 ) .* edge.ny;
FluxP(:, :, 2) = obj.gra .* fp( :, :, 6 ) .* edge.ny;

lambda = sqrt( max( obj.gra .* fm(:, :, 5), obj.gra .* fp(:, :, 5) ) );

FluxS(:, :, 1) = 0.5 .* ( FluxM(:, :, 1) + FluxP(:, :, 1) - ...
    lambda .* ( fp( :, :, 1 ) - fm( :, :, 1 ) ) );
FluxS(:, :, 2) = 0.5 .* ( FluxM(:, :, 2) + FluxP(:, :, 2) - ...
    lambda .* ( fp( :, :, 2 ) - fm( :, :, 2 ) ) );

frhs3d = edge.matEvaluateStrongFormEdgeRHS( FluxM, FluxS );

end

