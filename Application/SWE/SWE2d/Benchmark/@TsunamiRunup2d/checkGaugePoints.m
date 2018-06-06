function checkGaugePoints( obj )

tsunamiPos = makeNdgPostProcessFromNdgPhys( obj );
[ path, ~, ~ ] = fileparts( mfilename('fullpath') );
file = [path, '/mesh/output_ch5-7-9.xls'];
data = xlsread(file);
data(:, [2,3,4]) = data(:, [2,3,4])./100;

xg = [4.521, 4.521, 4.521];
yg = [1.196, 1.696, 2.196];
gaugePhys = tsunamiPos.interpolatePhysFieldToGaugePoint( obj.fphys, xg, yg, yg );
bot = gaugePhys(:, 4);

gaugeResult = tsunamiPos.interpolateOutputResultToGaugePoint(xg, yg, yg);

markersize = 6;
linewidth = 3;
fontsize = 14;

time = tsunamiPos.time{1};
% Gauge point #ch5
figure('color', 'w', 'Position', [100, 170, 625, 635]);
subplot(3,1,1); hold on; grid on; box on;
plot( data(:, 1), data(:, 2), 'ro', 'MarkerSize', markersize);
plot( time, gaugeResult( :, 1, 1 ) + bot(1), 'b', 'LineWidth', linewidth);
legend({'Measured', 'Numerical'}, 'box', 'off', 'Location', 'NorthWest', ...
    'FontSize', fontsize);
xlim( [ time(1), time(end) ] );

% Gauge point #ch7
subplot(3,1,2); hold on; grid on; box on;
plot(data(:, 1), data(:, 3), 'ro', 'MarkerSize', markersize);
plot( time, gaugeResult( :, 1, 2 ) + bot(2), 'b', 'LineWidth', linewidth);
xlim([0, obj.ftime]);

% Gauge point #ch9
subplot(3,1,3); hold on; grid on; box on;
plot(data(:, 1), data(:, 4), 'ro', 'MarkerSize', markersize);
plot( time, gaugeResult( :, 1, 3 ) + bot(3), 'b', 'LineWidth', linewidth);
xlim([0, obj.ftime]);
end

