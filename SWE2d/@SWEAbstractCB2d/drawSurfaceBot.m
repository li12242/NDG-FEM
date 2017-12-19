function drawSurfaceBot( obj )

for m = 1:obj.Nmesh
    mesh = obj.meshUnion(m);
    drawBottom( mesh, obj.fphys{m}(:,:,4) );
    eta = obj.fphys{m}(:,:,4) + obj.fphys{m}(:,:,1);
%     eta( obj.fphys{m}(:,:,1) < obj.hmin ) = nan;
    drawSurface( mesh, eta );
    %mesh.draw( eta );
    %set( mesh.figureHandle, 'FaceAlpha', 0.8 );
end
grid on;
view([20, 40])
end

function handle = drawSurface( mesh, eta )
% check the figure is created and not removed.
isFigureExit = ~isempty( mesh.figureHandle ) && isvalid( mesh.figureHandle ) ;
if ~isFigureExit
    handle = drawNew3dBottom( mesh, eta(:), [0, 0.2857, 0.8571] );
    set( handle, 'FaceAlpha', 0.8 );
end

end

function handle = drawBottom( mesh, bot )

% check the figure is created and not removed.
isFigureExit = ~isempty( mesh.figureHandle ) && isvalid( mesh.figureHandle ) ;
if ~isFigureExit
    figure;
    handle = drawNew3dBottom( mesh, bot(:), [.8, .8, .8] );
end

end% func

function handle = drawNew3dBottom( mesh, zvar, color )
EToV = ones(mesh.K, 1)*mesh.cell.Fmask(:)';
EToV = EToV + ( mesh.cell.Np*(0:mesh.K-1) )' * ones(1, mesh.cell.TNfp);
handle = patch(...
    'Vertices', [mesh.x(:), mesh.y(:), zvar(:)], ...
    'Faces', EToV, ...
    'FaceColor', color, ...
    'FaceVertexCData', zvar(:));
box on;
grid on;
end
