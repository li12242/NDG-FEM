function CheckGaugeResult( obj )

% tsunamiPos = makeNdgPostProcessFromNdgPhys( obj );
file = [ fileparts( mfilename('fullpath') ), ...
    '/mesh/output_ch5-7-9.xls'];
data = xlsread(file);
data(:, [2,3,4]) = data(:, [2,3,4])./100;

xg = [4.521, 4.521, 4.521];
yg = [1.196, 1.696, 2.196];
Ng = numel( xg );
analys = Analysis2d( obj, xg, yg );
% gaugePhys = tsunamiPos.interpolatePhysFieldToGaugePoint( obj.fphys, xg, yg, yg );
temp = analys.InterpGaugeResult( obj.fphys );
bot = temp(:, 4);
% gaugeResult = tsunamiPos.interpolateOutputResultToGaugePoint(xg, yg, yg);
Nout = obj.outputFile.outputStep;
gaugeResult = zeros( Nout, Ng );

for t = 1:Nout
    temp = analys.InterpOutputResult2GaugePoint( t );
    gaugeResult(t, :) = temp( :, 1 );
end

markersize = 6;
linewidth = 3;
fontsize = 14;

time = ncread( obj.outputFile.filename, 'time' );
% Gauge point #ch5
figure('color', 'w', 'Position', [100, 170, 625, 635]);
subplot(3,1,1); hold on; grid on; box on;
plot( data(:, 1), data(:, 2), 'ro', 'MarkerSize', markersize);
plot( time, gaugeResult( :, 1 ) + bot(1), 'b', 'LineWidth', linewidth);
legend({'Measured', 'Numerical'}, 'box', 'off', 'Location', 'NorthWest', ...
    'FontSize', fontsize);
xlim( [ time(1), time(end) ] );

% Gauge point #ch7
subplot(3,1,2); hold on; grid on; box on;
plot(data(:, 1), data(:, 3), 'ro', 'MarkerSize', markersize);
plot( time, gaugeResult( :, 2 ) + bot(2), 'b', 'LineWidth', linewidth);
xlim([0, obj.ftime]);

% Gauge point #ch9
subplot(3,1,3); hold on; grid on; box on;
plot(data(:, 1), data(:, 4), 'ro', 'MarkerSize', markersize);
plot( time, gaugeResult( :, 3 ) + bot(3), 'b', 'LineWidth', linewidth);
xlim([0, obj.ftime]);
end

