function drawMaxGaugeDepth( obj )

elv_max = [84.2, 49.1, 54.0, 40.2, 34.9, 27.4, 21.5, 16.1, 12.9];
bot_elv = [43.9, 34.5, 30.0, 27.4, 23.1, 19.1, 11.4, 09.3, 07.5];
dep_max = elv_max - bot_elv; % the maximum depth of each gauge points
xd = [5550, 11900, 13000, 4947.46, 5717.30, 6775.14,...
    7128.20, 8585.3, 9674.97, 10939.15, 11724.37, 12723.70...
    4913.11 5159.75 5790.63 5886.54 6763.05 6929.97 7326.02 7441.01...
    8735.94 8628.6 9761.13 9800 10957 11156.99 11689.05 11626.05 12333.72];

yd = [4400, 3250, 2700, 4289.71, 4407.61, 3869.23,...
    3162.00, 3443.08, 3085.89, 3044.78, 2810.41, 2485.08...
    4244.01 4369.62 4177.76 4503.97 3429.6 3591.87 2948.78 3232.12 3264.61...
    3604.63 3480.36 2414.79 2651.94 3800.72 2592.36 3406.8 2269.74];

malpos = makeNdgPostProcessFromNdgPhys( obj );
[ result ] = malpos.interpolateOutputResultToGaugePoint( xd, yd, xd );
[ gaugeValue ] = malpos.interpolatePhysFieldToGaugePoint( obj.fphys, xd, yd, xd );
bot = gaugeValue(:, 4)';
ind = 4:12; Ng = numel(ind);
dep_num_max = zeros(1, Ng);
for i = 1:9
    temp = result(ind(i), 1, :);
    dep_num_max(i) = max( temp(:) );
end

markersize = 8;
linewidth = 1.5;
fontsize = 12;

% draw the maximum depth at gauge points
figure('color', 'w');  hold on; box on;
plot(ind - 3, dep_max + bot_elv, 'ks-', ....
    'MarkerSize', markersize, 'LineWidth', linewidth);
plot(ind - 3, dep_num_max + bot_elv, 'ro--', ...
    'MarkerSize', markersize, 'LineWidth', linewidth);
legend({'Measured', 'Numerical'}, ...
    'box', 'off', 'Interpreter', 'Latex', 'FontSize', fontsize);

% draw the water depth at 13~29
elev = [79.15, 87.2, 54.9, 64.7, 51.1, 43.75, 44.35, 38.6, ...
    31.9, 40.75, 24.15, 24.9, 17.25, 20.7, 18.6, 17.25, 14]';
ind = 13:29; Ng = numel(ind);
dep_num_max = zeros(1, Ng);
for i = 1:9
    temp = result(ind(i), 1, :);
    dep_num_max(i) = max( temp(:) );
end
figure('color', 'w');  hold on; box on;
plot(ind - 12, elev, 'ks-', ....
    'MarkerSize', markersize, 'LineWidth', linewidth);
plot(ind - 12, dep_num_max + bot(ind), 'ro--', ...
    'MarkerSize', markersize, 'LineWidth', linewidth);
legend({'Measured', 'Numerical'}, ...
    'box', 'off', 'Interpreter', 'Latex', 'FontSize', fontsize);

% draw the time divergence
dep = zeros(malpos.Nt, 3);
for i = 1:3
    temp = result(i, 1, :);
    dep(:, i) = temp(:);
end
outputTime = ncread( malpos.outputFile{1}, 'time' );
time = zeros(3,1);
ind_A = find( dep(:, 1) > obj.hmin, 1 ); time(1) = outputTime(ind_A);
ind_B = find( dep(:, 2) > obj.hmin, 1 ); time(2) = outputTime(ind_B);
ind_C = find( dep(:, 3) > obj.hmin, 1 ); time(3) = outputTime(ind_C);

figure('color', 'w'); hold on; box on;
plot(1:2, [1140, 1320], ...
    'ks', 'MarkerFaceColor', 'k', 'MarkerSize', markersize);
plot(1:2, diff(time), ...
    'ro', 'MarkerFaceColor', 'r', 'MarkerSize', markersize);
legend({'Measured', 'Numerical'}, ...
    'box', 'off', 'Interpreter', 'Latex', 'FontSize', fontsize);

end% function