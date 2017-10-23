function draw(obj, varargin)

if( nargin == 1 )
    is2dFigure = true;
elseif( nargin == 2 )
    zvar = varargin{1};
    is2dFigure = false;
end
% check the figure is created and not removed.
isFigureExit = ~isempty( obj.figureHandle ) && isvalid( obj.figureHandle ) ;

if isFigureExit
    if ~is2dFigure % update 3d figure
        updateFigure(obj.figureHandle, obj.x(:), obj.y(:), zvar(:));
    end
else
    if is2dFigure
        obj.figureHandle = drawNew2dFigure(obj);
    else
        obj.figureHandle = drawNew3dFigure(obj, zvar(:) );
    end
end

end% func

function handle = drawNew2dFigure(obj)
figure;
handle = patch(...
    'Vertices', [obj.vx(:), obj.vy(:)], ...
    'Faces', obj.EToV', ...
    'FaceColor', [0.8, 0.9, 1]);
box on;
grid on;
end

function handle = drawNew3dFigure(obj, zvar)
EToV = ones(obj.K, 1)*obj.cell.Fmask(:)';
EToV = EToV + ( obj.cell.Np*(0:obj.K-1) )' * ones(1, obj.cell.TNfp);
handle = patch(...
    'Vertices', [obj.x(:), obj.y(:), zvar(:)], ...
    'Faces', EToV, ...
    'FaceColor', 'interp', ...
    'FaceVertexCData', zvar(:));
box on;
grid on;
end

function updateFigure(handle, xvar, yvar, zvar)
set(handle, 'Vertices', [xvar(:), yvar(:), zvar(:)],...
    'FaceVertexCData', zvar(:));
end