function drawSurfaceBot( obj )

for m = 1:obj.Nmesh
    mesh = obj.meshUnion(m);
    drawBottom( mesh, obj.fphys{m}(:,:,4) );
    eta = obj.fphys{m}(:,:,4) + obj.fphys{m}(:,:,1);
    eta( obj.fphys{m}(:,:,1) < obj.hmin ) = nan;
    mesh.draw( eta );
end
grid on;
view([20, 40])
end

function drawBottom( mesh, bot )

% check the figure is created and not removed.
isFigureExit = ~isempty( mesh.figureHandle ) && isvalid( mesh.figureHandle ) ;
if ~isFigureExit
    drawNew3dBottom( mesh, bot(:) );
end

end% func

function handle = drawNew3dBottom( mesh, zvar)
EToV = ones(mesh.K, 1)*mesh.cell.Fmask(:)';
EToV = EToV + ( mesh.cell.Np*(0:mesh.K-1) )' * ones(1, mesh.cell.TNfp);
handle = patch(...
    'Vertices', [mesh.x(:), mesh.y(:), zvar(:)], ...
    'Faces', EToV, ...
    'FaceColor', [.8, .8, .8], ...
    'FaceVertexCData', zvar(:));
box on;
grid on;
end
