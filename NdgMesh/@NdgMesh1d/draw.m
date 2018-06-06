function draw( obj, zvar )
% check the figure is created and not removed.
isFigureExit = ~isempty( obj.figureHandle ) ...
    && isvalid( obj.figureHandle ) ;

if isFigureExit
    %updateFigure(obj.figureHandle, obj.x(:), zvar(:));
    set( obj.figureHandle, 'YData', zvar(:) );
else
    Np = obj.cell.Np;
    K = obj.K;
    list = 1:Np:(K*Np);
    g = graph();
    for n = 1:Np-1
        g = addedge(g, list+n-1, list+n);
    end
    
    obj.figureHandle = plot(g, ...
        'XData', obj.x(:), 'YData', zvar(:), ...
        'LineWidth', 1, ...
        'Marker', 'o', ...
        'NodeColor','k', ...
        'EdgeColor', 'k', ...
        'MarkerSize', 2, ...
        'NodeLabel', {});
    box on;
    grid on;
end

end