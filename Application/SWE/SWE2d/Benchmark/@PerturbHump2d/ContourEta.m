function ContourEta( obj )
%CONTOURETA Summary of this function goes here
%   Detailed explanation goes here

time = [0.12, 0.24, 0.36, 0.48, 0.6];
Nstep = numel( time );

Ng = 300; tol = 0;
xg = linspace( 0+tol, 2-tol, Ng );
yg = linspace( 0+tol, 1-tol, Ng );
[ x, y ] = meshgrid( xg, yg );
analysis = Analysis2d( obj, x(:), y(:) );
data = analysis.InterpGaugeResult( obj.fphys );
%v2 = [linspace( 0.999, 1.0007, 30 )];
%v1 = [linspace( 0.999, 1.0003, 30 )];
for n = 1:Nstep
    step = find( obj.outputFile.outputTime >= time(n), 1 );

    figure;
    flux = analysis.InterpOutputResult2GaugePoint( step );
    contourf( x, y, real( reshape( flux(:,1) + data(:,4), Ng, Ng) ), 20 );

    if n < 3
        set(gca, 'CLim', [0.998, 1.008])
    elseif n == 3
        set(gca, 'CLim', [0.996, 1.008])
    else
        set(gca, 'CLim', [0.996, 1.006])
    end
    colorbar;
    colormap jet;
    
    xlabel('$x$(m)', 'FontSize', 16, 'Interpreter', 'Latex');
    ylabel('$y$(m)', 'FontSize', 16, 'Interpreter', 'Latex');
end

end

