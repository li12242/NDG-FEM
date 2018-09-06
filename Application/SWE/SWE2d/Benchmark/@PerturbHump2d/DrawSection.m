function DrawSection( obj )

time = 0.2;

Ng = 400; tol = 0;
xg = linspace( 0+tol, 2-tol, Ng );
yg = zeros( size(xg) );
analysis = Analysis2d( obj, xg(:), yg(:) );
data = analysis.InterpGaugeResult( obj.fphys );
%v2 = [linspace( 0.999, 1.0007, 30 )];
%v1 = [linspace( 0.999, 1.0003, 30 )];
step = find( obj.outputFile.outputTime >= time, 1 );
figure;
flux = analysis.InterpOutputResult2GaugePoint( step );

subplot(1, 2, 1);
subplot(1, 2, 2);

xlabel('$x$(m)', 'FontSize', 16, 'Interpreter', 'Latex');
ylabel('$y$(m)', 'FontSize', 16, 'Interpreter', 'Latex');

end