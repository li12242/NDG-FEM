function CoutourFlux( obj )
Ng = 50; tol = 0;
xg = linspace( 0+tol, 1-tol, Ng );
yg = linspace( 0+tol, 1-tol, Ng );
[ x, y ] = meshgrid( xg, yg );
analysis = Analysis2d( obj, x(:), y(:) );
flux = analysis.InterpOutputResult2GaugePoint( obj.outputFile.outputStep );
contourf( x, y, reshape( real( flux(:,2) ), Ng, Ng), 20 );
colorbar;
colormap jet;
xlabel('$x$(m)', 'FontSize', 16, 'Interpreter', 'Latex');
ylabel('$y$(m)', 'FontSize', 16, 'Interpreter', 'Latex');
end