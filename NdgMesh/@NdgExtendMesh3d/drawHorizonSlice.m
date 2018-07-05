function drawHorizonSlice( obj, nodeField3d )

Nhp = obj.mesh2d.cell.Np;
Nzp = obj.cell.Nz + 1;

isFigureExit = ~isempty( obj.figure_handle );
if isFigureExit
    isFigureExit = isvalid( obj.figure_handle{1} );
else
    obj.figure_handle = cell( Nzp, 1);
end

if isFigureExit
    for n = 1 : Nzp
        ind = (1 : Nhp) + Nhp * ( n - 1 );
        Update_Figure(obj.figure_handle{n}, ...
            obj.x(ind, :), obj.y(ind, :), obj.z(ind, :), nodeField3d(ind, :));
    end
else
    for n = 1 : Nzp
        ind = (1 : Nhp) + Nhp * ( n - 1 );
        obj.figure_handle{n} = drawNew3dFigure(obj, ...
            obj.x(ind, :), obj.y(ind, :), obj.z(ind, :), nodeField3d(ind, :) );
    end
end

end% func

function handle = drawNew3dFigure(obj, x, y, z, cvar)
mesh2d = obj.mesh2d;
EToV = ones(obj.K, 1)*mesh2d.cell.Fmask(:)';
EToV = EToV + ( mesh2d.cell.Np*( 0 : obj.K-1 ) )' * ones(1, mesh2d.cell.TNfp);
handle = patch(...
    'Vertices', [x(:), y(:), z(:)], ...
    'Faces', EToV, ...
    'FaceColor', 'interp', ...
    'FaceVertexCData', cvar(:));
box on;
grid on;
end

function Update_Figure(handle, xvar, yvar, zvar, cvar)
set(handle, 'Vertices', [xvar(:), yvar(:), zvar(:)],...
    'FaceVertexCData', cvar(:));
end