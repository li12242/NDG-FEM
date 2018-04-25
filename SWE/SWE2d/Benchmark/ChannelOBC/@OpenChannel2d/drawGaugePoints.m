function drawGaugePoints( obj, xg )

isDrawExact = 1; % draw exact surface line
Ng = numel(xg);
%xg = linspace( 0, obj.ChLength, Ng );
yg = ones( size( xg ) ) * obj.ChWidth / 2;

pos = makeNdgPostProcessFromNdgPhys( obj );
[ gaugeValue ] = pos.interpolateOutputResultToGaugePoint( xg, yg, yg );
[ time ] = pos.time{1};

figure;
for n = 1:Ng
    subplot(Ng, 1, n); 
    hold on; grid on; box on;
    if isDrawExact
        [ hext, ~ ] = obj.setExactSolution( xg(n), time );
        plot( time, hext(:) - obj.H, 'r:', 'LineWidth', 2 );
    end
    temp = gaugeValue(:, 1, n);
    plot( time, temp(:) - obj.H, 'b-', 'LineWidth', 1.5 );
    %ylim([-.25, .25]);
%     xlim([time(end) - 2000, time(end)]);
    [ pxx, f ] = plomb( temp(1:end-1), time(1:end-1), [], 10, 'power' );
    [ pk, f0 ] = findpeaks( pxx, f, 'MinPeakHeight', 0.01);
    fprintf(' No.%d,  eta = %e, T = %e\n', n, pk, 1./f0);
end

end
