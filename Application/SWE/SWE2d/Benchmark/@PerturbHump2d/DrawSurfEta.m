function DrawSurfEta( obj )
%DRAWSURFETA Summary of this function goes here
%   Detailed explanation goes here

time = [0.12, 0.24, 0.36, 0.48, 0.6];
Nstep = numel( time );

% Ng = 100; tol = 0;
% xg = linspace( 0+tol, 2-tol, Ng );
% yg = linspace( 0+tol, 1-tol, Ng );
% [ x, y ] = meshgrid( xg, yg );
% analysis = Analysis2d( obj, x(:), y(:) );
% data = analysis.InterpGaugeResult( obj.fphys );
for n = 1:Nstep
    step = find( obj.outputFile.outputTime >= time(n), 1 );

    figure;
    fphys = obj.outputFile.readOutputResult( step );
    visual = Visual2d( obj.meshUnion );
    visual.drawResult( fphys(:, :, 1) + obj.fphys{1}(:, :, 4) );
    visual.drawHandle.EdgeColor = 'none';
    set(gca, 'CLim', [0.996, 1.006]);
    colormap jet;
    xlabel('$x$(m)', 'FontSize', 16, 'Interpreter', 'Latex');
    ylabel('$y$(m)', 'FontSize', 16, 'Interpreter', 'Latex');
end

end

