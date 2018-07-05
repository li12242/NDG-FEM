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
        obj, mesh2d.BoundaryEdge, fphys2d, obj.fext2d{m});
    
    % extend 3d field
    fphys3d = Evaluate3d_Auxiliary( obj, mesh3d, fphys2d, fphys3d );
    
    % evaluate 3d velocity volume integral
    obj.frhs3d{m} = Evaluate3d_VolumeKernel( obj, obj.gra, mesh3d, fphys3d{m} );
    
    % evaluate 3d velocity horizontal surface integral
    obj.frhs3d{m} = obj.frhs3d{m} + Evaluate3d_SideSurfaceKernel( ...
        obj.gra, mesh3d.InnerEdge, fphys3d );
    
    obj.frhs3d{m} = obj.frhs3d{m} + Evaluate3d_HorizontalBoundaryKernel( ...
        obj, mesh3d.BoundaryEdge, fphys3d, obj.fext3d{m} );
    
    % evaluate 3d velocity field surface integral
    obj.frhs3d{m} = obj.frhs3d{m} + Evaluate3d_SurfaceBoundaryKernel( ...
        obj, mesh3d.SurfaceBoundaryEdge, fphys3d );
    
    % evaluate 3d velocity field bottom integral
    obj.frhs3d{m} = obj.frhs3d{m} + Evaluate3d_BottomSurfaceKernel( ...
        obj, mesh3d.BottomEdge, fphys3d);
    
    obj.frhs3d{m} = obj.frhs3d{m} + Evaluate3d_BottomBoundaryKernel( ...
        obj, mesh3d.BottomBoundaryEdge, fphys3d );
end

end

function fphys2d = EvaluateHorizon2d_Kernel( mesh3d, fphys2d, fphys3d )
% evaluate 2d averaged velocity
U = mesh3d.VerticalColumnIntegralField( fphys3d(:, :, 1) );
V = mesh3d.VerticalColumnIntegralField( fphys3d(:, :, 2) );

%fphys2d(:, :, 4) = fphys2d(:, :, 1) - fphys2d(:, :, 5);
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

function fphys3d = Evaluate3d_Auxiliary( obj, mesh3d, fphys2d, fphys3d )

fphys3d{1}(:, :, 6) = mesh3d.Extend2dField( fphys2d{1}(:, :, 4) );
fphys3d{1}(:, :, 7) = mesh3d.Extend2dField( fphys2d{1}(:, :, 1) );
% fphys3d(:, :, 5) = mesh3d.Extend2dField( fphys2d(:, :, 4) );
% fphys3d(:, :, 6) = mesh3d.Extend2dField( fphys2d(:, :, 1) );

% Hmiu = miu ./ fphys3d{1}(:, :, 6) .^ 2;
Hmiu = sqrt( obj.miu0 );
% Hmiu = miu;

% evaluate 3d auxiliary variable
edge3d = mesh3d.BottomEdge;
[ fm, fp ] = edge3d.matEvaluateSurfValue( fphys3d );
FluxM_1(:, :, 1) = fm(:, :, 1); % normal vector is - 1
FluxM_1(:, :, 2) = fm(:, :, 2);
FluxP_1(:, :, 1) = fp(:, :, 1);
FluxP_1(:, :, 2) = fp(:, :, 2);
% FluxS_1 = 0.5 *( FluxM_1 + FluxP_1 );
% FluxS_1 = FluxP_1;
% fphys3d{1}(:, :, 4:5) = edge3d.matEvaluateStrongFormEdgeRHS( FluxM_1, FluxP_1, FluxS_1 );
% fphys3d{1}(:, :, 4:5) = edge3d.matEvaluateStrongFormEdgeAlterRHS( FluxM_1, FluxP_1 );
fphys3d{1}(:, :, 4:5) = edge3d.matEvaluateStrongFormEdgeCentralRHS( FluxM_1, FluxP_1 );

edge3d = mesh3d.BottomBoundaryEdge;
[ fm, fp ] = edge3d.matEvaluateSurfValue( fphys3d );
FluxM(:, :, 1) = fm(:, :, 1);
FluxM(:, :, 2) = fm(:, :, 2);
FluxS(:, :, 1) = fp(:, :, 1);
FluxS(:, :, 2) = fp(:, :, 2);
fphys3d{1}(:, :, 4:5) = fphys3d{1}(:, :, 4:5) ...
    + edge3d.matEvaluateStrongFormEdgeRHS( FluxM, FluxS );

fphys3d{1}(:, :, 4) = - Hmiu .* ( fphys3d{1}(:, :, 4) + ...
    mesh3d.tz .* ( mesh3d.cell.Dt * fphys3d{1}(:,:,1) ) );

fphys3d{1}(:, :, 5) = - Hmiu .* ( fphys3d{1}(:, :, 5) + ...
    mesh3d.tz .* ( mesh3d.cell.Dt * fphys3d{1}(:,:,2) ) );

% fphys3d{1}(:, :, 4) = - Hmiu .* mesh3d.tz .* ( mesh3d.cell.Dt * fphys3d{1}(:,:,1) );
% fphys3d{1}(:, :, 5) = - Hmiu .* mesh3d.tz .* ( mesh3d.cell.Dt * fphys3d{1}(:,:,2) );
end

function frhs3d = Evaluate3d_VolumeKernel( obj, gra, mesh3d, fphys3d )

phi = gra * fphys3d(:, :, 7);

dFdr = mesh3d.cell.Dr * phi;
dFds = mesh3d.cell.Ds * phi;
% dFdt = mesh3d.cell.Dt * phi;
Hmiu = sqrt( obj.miu0 );

frhs3d(:, :, 1) = - ( mesh3d.rx .* dFdr + mesh3d.sx .* dFds + ...
    Hmiu .* (mesh3d.tz .* ( mesh3d.cell.Dt * fphys3d(:, :, 4) ) ) );

frhs3d(:, :, 2) = - ( mesh3d.ry .* dFdr + mesh3d.sy .* dFds + ...
    Hmiu .* (mesh3d.tz .* ( mesh3d.cell.Dt * fphys3d(:, :, 5) ) ) );
end

function frhs3d = Evaluate3d_SideSurfaceKernel( gra, edge, fphys3d )

[ fm, fp ] = edge.matEvaluateSurfValue( fphys3d );
FluxM(:, :, 1) = gra .* fm( :, :, 7 ) .* edge.nx;
FluxP(:, :, 1) = gra .* fp( :, :, 7 ) .* edge.nx;
FluxM(:, :, 2) = gra .* fm( :, :, 7 ) .* edge.ny;
FluxP(:, :, 2) = gra .* fp( :, :, 7 ) .* edge.ny;

lambda = max( sqrt( gra .* fm(:, :, 6) ), sqrt( gra .* fp(:, :, 6) ) );

FluxS(:, :, 1) = 0.5 .* ( FluxM(:, :, 1) + FluxP(:, :, 1) - ...
    lambda .* ( fp( :, :, 1 ) - fm( :, :, 1 ) ) );
FluxS(:, :, 2) = 0.5 .* ( FluxM(:, :, 2) + FluxP(:, :, 2) - ...
    lambda .* ( fp( :, :, 2 ) - fm( :, :, 2 ) ) );

frhs3d = edge.matEvaluateStrongFromEdgeRHS( FluxM, FluxP, FluxS );
end

function frhs3d = Evaluate3d_SurfaceBoundaryKernel( obj, edge, fphys3d )

Hmiu = sqrt( obj.miu0 );

[ fm, ~ ] = edge.matEvaluateSurfValue( fphys3d );
FluxM(:, :, 1) = Hmiu .* fm( :, :, 4 ) .* edge.nz;
FluxM(:, :, 2) = Hmiu .* fm( :, :, 5 ) .* edge.nz;
FluxS(:, :, 1) = - Hmiu .* fm( :, :, 4 ) .* edge.nz;
FluxS(:, :, 2) = - Hmiu .* fm( :, :, 5 ) .* edge.nz;

frhs3d = edge.matEvaluateStrongFormEdgeRHS( FluxM, FluxS );
end

function frhs3d = Evaluate3d_BottomSurfaceKernel( ...
    obj, edge, fphys3d )

Hmiu = sqrt( obj.miu0 );
tau = 1e-3;

[ fm, fp ] = edge.matEvaluateSurfValue( fphys3d );
FluxM(:, :, 1) = Hmiu .* ( fm( :, :, 4 ) .* edge.nz + tau * fm(:, :, 1) );
FluxM(:, :, 2) = Hmiu .* ( fm( :, :, 5 ) .* edge.nz + tau * fm(:, :, 2) );
FluxP(:, :, 1) = Hmiu .* ( fp( :, :, 4 ) .* edge.nz + tau * fp(:, :, 1) );
FluxP(:, :, 2) = Hmiu .* ( fp( :, :, 5 ) .* edge.nz + tau * fp(:, :, 2) );

% FluxS(:, :, 1) = 0.5 .* ( FluxM(:, :, 1) + FluxP(:, :, 1)  );
% FluxS(:, :, 2) = 0.5 .* ( FluxM(:, :, 2) + FluxP(:, :, 2)  );
% FluxS(:, :, 1) = FluxP(:, :, 1);
% FluxS(:, :, 2) = FluxP(:, :, 2);

% frhs3d = edge.matEvaluateStrongFormEdgeRHS( FluxM, FluxP, FluxS );
frhs3d = edge.matEvaluateStrongFormEdgeCentralRHS( FluxM, FluxP );
end

function frhs3d = Evaluate3d_BottomBoundaryKernel( obj, edge, fphys3d )
Hmiu = sqrt( obj.miu0 );

[ fm, fp ] = edge.matEvaluateSurfValue( fphys3d );
FluxM(:, :, 1) = Hmiu .* fm( :, :, 4 ) .* edge.nz;
FluxM(:, :, 2) = Hmiu .* fm( :, :, 5 ) .* edge.nz;

FluxS(:, :, 1) =  - obj.K .* fp( :, :, 1 ) ./ fp( :, :, 6 ) .* edge.nz;
FluxS(:, :, 2) =  - obj.K .* fp( :, :, 2 ) ./ fp( :, :, 6 ) .* edge.nz;

frhs3d = edge.matEvaluateStrongFormEdgeRHS( FluxM, FluxS );
end
