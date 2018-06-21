function matEvaluateRHS( obj, fphys2d, fphys3d )
%MATEVALUATERHS Summary of this function goes here
%   Detailed explanation goes here

for m = 1:obj.Nmesh
    mesh2d = obj.mesh2d(m);
    mesh3d = obj.mesh3d(m);
    
    % evaluate 2d PCE volume integral term
    fphys2d{m} = EvaluateHorizon2d_Kernel( mesh3d, ...
        fphys2d{m}, fphys3d{m} );
    
    obj.frhs2d{m} = EvaluatePCE2d_VolumeKernel( mesh2d, fphys2d{m} );

    % evaluate 2d PCE surface integral term
    obj.frhs2d{m} = obj.frhs2d{m} + EvaluatePCE2d_SurfaceKernel( ...
        obj.gra, mesh2d.InnerEdge, fphys2d );
    
    obj.frhs2d{m} = obj.frhs2d{m} + EvaluatePCE2d_BoundaryKernel( ...
        mesh2d.BoundaryEdge, fphys2d);
    
    % extend 3d field
    fphys3d{m} = Evaluate3d_Auxiliary( obj.miu, mesh3d, fphys2d{m}, fphys3d{m} );
    
    % evaluate 3d velocity volume integral
    obj.frhs3d{m} = Evaluate3d_VolumeKernel( obj.gra, mesh3d, fphys3d{m} );
    
    % evaluate 3d velocity surface integral
    obj.frhs3d{m} = obj.frhs3d{m} + Evaluate3d_SideSurfaceKernel( ...
        obj.gra, mesh3d.InnerEdge, fphys3d );
    
    obj.frhs3d{m} = obj.frhs3d{m} + Evaluate3d_BoundarySurfaceKernel( ...
        obj.gra, mesh3d.BoundaryEdge, fphys3d );
    
    % evaluate 3d velocity field bottom surface integral
%     obj.frhs3d{m} = obj.frhs3d{m} + Evaluate3d_BottomSurfaceKernel( ...
%         obj.gra, obj.K, mesh3d.BottomEdge, fphys3d);
    
end

end

function frhs2d = EvaluatePCE2d_BoundaryKernel( edge, fphys2d )
[ fm, fp ] = edge.matEvaluateSurfValue( fphys2d );

un =  fm( :, :, 2 ) .* edge.nx + fm( :, :, 3).* edge.ny;
vn = -fm( :, :, 2 ) .* edge.ny + fm( :, :, 3).* edge.nx;

fp(:, :, 2) = - un .* edge.nx - vn .* edge.ny;
fp(:, :, 3) = - un .* edge.ny + vn .* edge.nx;

FluxM = fm(:, :, 2) .* edge.nx + fm(:, :, 3) .* edge.ny;
FluxP = fp(:, :, 2) .* edge.nx + fp(:, :, 3) .* edge.ny;
FluxS = 0.5 * ( FluxM + FluxP );

frhs2d = edge.matEvaluateStrongFromEdgeRHS( FluxM, FluxS );
end

function fphys2d = EvaluateHorizon2d_Kernel( mesh3d, fphys2d, fphys3d )
% evaluate 2d averaged velocity
U = mesh3d.VerticalIntegralField( fphys3d(:, :, 1) );
V = mesh3d.VerticalIntegralField( fphys3d(:, :, 2) );

fphys2d(:, :, 4) = fphys2d(:, :, 1) - fphys2d(:, :, 5);
fphys2d(:, :, 2) = fphys2d(:, :, 4) .* U;
fphys2d(:, :, 3) = fphys2d(:, :, 4) .* V;
end

function frhs2d = EvaluatePCE2d_VolumeKernel( mesh2d, fphys2d )

% evaluated 2d depth rhs
frhs2d(:, :, 1) = -( ...
    mesh2d.rx .* ( mesh2d.cell.Dr * fphys2d(:, :, 2) ) + ...
    mesh2d.sx .* ( mesh2d.cell.Ds * fphys2d(:, :, 2) ) + ...
    mesh2d.ry .* ( mesh2d.cell.Dr * fphys2d(:, :, 3) ) + ...
    mesh2d.sy .* ( mesh2d.cell.Ds * fphys2d(:, :, 3) ) );

end

function frhs2d = EvaluatePCE2d_SurfaceKernel( gra, edge, fphys2d )

[ fm, fp ] = edge.matEvaluateSurfValue( fphys2d );
lambda = max( sqrt( gra .* fm(:, :, 4) ), sqrt( gra .* fp(:, :, 4) ) );

FluxM = fm(:, :, 2) .* edge.nx + fm(:, :, 3) .* edge.ny;
FluxP = fp(:, :, 2) .* edge.nx + fp(:, :, 3) .* edge.ny;
FluxS = 0.5 * ( FluxM + FluxP - lambda .* ( fp(:, :, 1) - fm(:, :, 1)  ) );

frhs2d = edge.matEvaluateStrongFromEdgeRHS( FluxM, FluxP, FluxS );
end

function fphys3d = Evaluate3d_Auxiliary( miu, mesh3d, fphys2d, fphys3d )

fphys3d(:, :, 5) = mesh3d.Extend2dField( fphys2d(:, :, 4) );
Hmiu = miu ./ fphys3d(:, :, 5) .^ 2;
fphys3d(:, :, 6) = mesh3d.Extend2dField( fphys2d(:, :, 1) );

% evaluate 3d auxiliary variable
fphys3d(:, :, 3) = - Hmiu .* ( mesh3d.tz .* ( mesh3d.cell.Dt * fphys3d(:,:,1) ) );
fphys3d(:, :, 4) = - Hmiu .* ( mesh3d.tz .* ( mesh3d.cell.Dt * fphys3d(:,:,2) ) );

end

function frhs3d = Evaluate3d_VolumeKernel( gra, mesh3d, fphys3d )

phi = gra * fphys3d(:, :, 6);

dFdr = mesh3d.cell.Dr * phi;
dFds = mesh3d.cell.Ds * phi;
% dFdt = mesh3d.cell.Dt * phi;

frhs3d(:, :, 1) = - ( mesh3d.rx .* dFdr + mesh3d.sx .* dFds + ...
    mesh3d.tz .* ( mesh3d.cell.Dt * fphys3d(:, :, 3) ) );

frhs3d(:, :, 2) = - ( mesh3d.ry .* dFdr + mesh3d.sy .* dFds + ...
    mesh3d.tz .* ( mesh3d.cell.Dt * fphys3d(:, :, 4) ) );
end

function frhs3d = Evaluate3d_SideSurfaceKernel( gra, edge, fphys3d )

[ fm, fp ] = edge.matEvaluateSurfValue( fphys3d );
FluxM(:, :, 1) = gra .* fm( :, :, 6 ) .* edge.nx;
FluxP(:, :, 1) = gra .* fp( :, :, 6 ) .* edge.nx;
FluxM(:, :, 2) = gra .* fm( :, :, 6 ) .* edge.ny;
FluxP(:, :, 2) = gra .* fp( :, :, 6 ) .* edge.ny;

lambda = max( sqrt( gra .* fm(:, :, 5) ), sqrt( gra .* fp(:, :, 5) ) );

FluxS(:, :, 1) = 0.5 .* ( FluxM(:, :, 1) + FluxP(:, :, 1) - ...
    lambda .* ( fp( :, :, 1 ) - fm( :, :, 1 ) ) );
FluxS(:, :, 2) = 0.5 .* ( FluxM(:, :, 2) + FluxP(:, :, 2) - ...
    lambda .* ( fp( :, :, 2 ) - fm( :, :, 2 ) ) );

frhs3d = edge.matEvaluateStrongFromEdgeRHS( FluxM, FluxP, FluxS );
end

function frhs3d = Evaluate3d_BoundarySurfaceKernel( gra, edge, fphys3d )

[ fm, fp ] = edge.matEvaluateSurfValue( fphys3d );

un =  fm( :, :, 1 ) .* edge.nx + fm( :, :, 2).* edge.ny;
vn = -fm( :, :, 1 ) .* edge.ny + fm( :, :, 2).* edge.nx;

fp(:, :, 1) = - un .* edge.nx - vn .* edge.ny;
fp(:, :, 2) = - un .* edge.ny + vn .* edge.nx;

FluxM(:, :, 1) = gra .* fm( :, :, 6 ) .* edge.nx;
FluxP(:, :, 1) = gra .* fp( :, :, 6 ) .* edge.nx;
FluxM(:, :, 2) = gra .* fm( :, :, 6 ) .* edge.ny;
FluxP(:, :, 2) = gra .* fp( :, :, 6 ) .* edge.ny;

lambda = sqrt( gra .* fm(:, :, 5) );

FluxS(:, :, 1) = 0.5 .* ( FluxM(:, :, 1) + FluxP(:, :, 1) - ...
    lambda .* ( fp( :, :, 1 ) - fm( :, :, 1 ) ) );
FluxS(:, :, 2) = 0.5 .* ( FluxM(:, :, 2) + FluxP(:, :, 2) - ...
    lambda .* ( fp( :, :, 2 ) - fm( :, :, 2 ) ) );

frhs3d = edge.matEvaluateStrongFromEdgeRHS( FluxM, FluxS );

end

function frhs3d = Evaluate3d_BottomSurfaceKernel( ...
    gra, K, edge, fphys3d )

[ fm, fp ] = edge.matEvaluateSurfValue( fphys3d );
FluxM(:, :, 1) = fm( :, :, 3 ) .* edge.nz;
FluxM(:, :, 2) = fm( :, :, 4 ) .* edge.nz;
FluxP(:, :, 1) = fp( :, :, 3 ) .* edge.nz;
FluxP(:, :, 2) = fp( :, :, 4 ) .* edge.nz;

lambda = max( sqrt( gra .* fm(:, :, 5) ), sqrt( gra .* fp(:, :, 5)  ) );

FluxS(:, :, 1) = 0.5 .* ( FluxM(:, :, 1) + FluxP(:, :, 1) - ...
    lambda .* ( fp( :, :, 1 ) - fm( :, :, 1 ) ) );
FluxS(:, :, 2) = 0.5 .* ( FluxM(:, :, 2) + FluxP(:, :, 2) - ...
    lambda .* ( fp( :, :, 2 ) - fm( :, :, 2 ) ) );

ind = (edge.ftype == enumBoundaryCondition.BottomBoundary);
FluxS(:, ind, 1) = - K .* fp( :, ind, 1 );
FluxS(:, ind, 2) = - K .* fp( :, ind, 2 );

frhs3d = edge.matEvaluateStrongFormEdgeRHS( FluxM, FluxP, FluxS );
end
