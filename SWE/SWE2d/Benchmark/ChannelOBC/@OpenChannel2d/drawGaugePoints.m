function drawGaugePoints( obj, xg )

isDrawExact = 0; % draw exact surface line
Ng = 5;
%xg = linspace( 0, obj.ChLength, Ng );
yg = ones( size( xg ) ) * obj.ChWidth / 2;

pos = makeNdgPostProcessFromNdgPhys( obj );
[ gaugeValue ] = pos.interpolateOutputResultToGaugePoint( xg, yg, yg );
[ time ] = ncread( pos.outputFile{1}, 'time' );

for n = 1:Ng
    subplot(Ng, 1, n)
    if isDrawExact
        [ hext, ~ ] = obj.setExactSolution( xg(n), time );
        plot( time, hext(:) - obj.H, 'r--', 'LineWidth', 2.5 );
    end
    temp = gaugeValue(n, 1, :);
    temp = temp(:);
    hold on; grid on; box on;
    plot( time, temp(:) - obj.H, 'b--', 'LineWidth', 2 );
    ylim([-.45, .45]);
    xlim([time(end) - 2000, time(end)]);
    [ pxx, f ] = plomb( temp(1:end-1), time(1:end-1), [], 10, 'power' );
    [ pk, f0 ] = findpeaks( pxx, f, 'MinPeakHeight', 0.001);
    fprintf(' No.%d,  eta = %e, T = %e\n', n, pk, 1./f0);
end

end
