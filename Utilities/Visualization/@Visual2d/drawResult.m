function drawResult( obj, field )
    
    isFigureEmpty = isempty( obj.drawHandle );
    if isFigureEmpty
        drawNewFigure( obj, field );
        return;
    end

    if isvalid( obj.drawHandle )
        updateFigure( obj, field );
    else
        drawNewFigure( obj, field );
    end
end

function updateFigure( obj, field )
    set(obj.drawHandle, 'Vertices', [obj.mesh.x(:), obj.mesh.y(:), field(:)],...
        'FaceVertexCData', field(:));
end

function drawNewFigure(obj, field)
    obj.drawHandle = trisurf( obj.tri, obj.mesh.x(:), obj.mesh.y(:), field(:) );
    box on;
    grid on;
end
