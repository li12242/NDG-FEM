function draw( obj, varargin )

switch nargin
    case 1
        fphys = obj.fphys;
    case 2
        fphys = varargin{1};
end

for m = 1:obj.Nmesh
    mesh = obj.meshUnion(m);
    drawBottom( mesh, fphys{m}(:,:,4) );
    eta = fphys{m}(:,:,4) + fphys{m}(:,:,1);
%     eta( fphys{m}(:,:,1) < obj.hmin ) = nan;
    drawSurface( mesh, eta );
end
grid on;
view([20, 40])
% zlim([0, 6]);
colormap winter;
set( gca, 'CLim', [2.4, 5] );
xlabel('$x$ (m)', 'FontSize', 16, 'Interpreter', 'Latex');
ylabel('$y$ (m)', 'FontSize', 16, 'Interpreter', 'Latex');
zlabel('$\eta$ (m)', 'FontSize', 16, 'Interpreter', 'Latex');
end

function handle = drawSurface( mesh, eta )
% check the figure is created and not removed.
isFigureExit = ~isempty( mesh.figureHandle ) && isvalid( mesh.figureHandle ) ;
if ~isFigureExit
    handle = drawNew3dSurface( mesh, eta(:));%, [0, 0.2857, 0.8571] );
    set( handle, 'FaceAlpha', 0.8 );
end

end

function handle = drawBottom( mesh, bot )

% check the figure is created and not removed.
isFigureExit = ~isempty( mesh.figureHandle ) && isvalid( mesh.figureHandle ) ;
if ~isFigureExit
    figure;
    handle = drawNew3dBottom( mesh, bot(:), [.67, .67, .67] );
end

end% func

function handle = drawNew3dSurface( mesh, zvar )
EToV = ones(mesh.K, 1)*mesh.cell.Fmask(:)';
EToV = EToV + ( mesh.cell.Np*(0:mesh.K-1) )' * ones(1, mesh.cell.TNfp);
handle = patch(...
    'Vertices', [mesh.x(:), mesh.y(:), zvar(:)], ...
    'Faces', EToV, ...
    'FaceColor', 'interp', ...
    'EdgeColor', 'k', ...
    'FaceVertexCData', zvar(:));
box on;
grid on;
end

function handle = drawNew3dBottom( mesh, zvar, color )
EToV = ones(mesh.K, 1)*mesh.cell.Fmask(:)';
EToV = EToV + ( mesh.cell.Np*(0:mesh.K-1) )' * ones(1, mesh.cell.TNfp);
handle = patch(...
    'Vertices', [mesh.x(:), mesh.y(:), zvar(:)], ...
    'Faces', EToV, ...
    'FaceColor', color, ...
    'EdgeColor', [0.0588, 0.0196, 0.4], ...
    'FaceVertexCData', zvar(:));
box on;
grid on;
end
