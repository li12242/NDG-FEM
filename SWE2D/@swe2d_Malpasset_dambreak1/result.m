function result( obj )
%RESULT Summary of this function goes here
%   Detailed explanation goes here

% the maximum water elevation at gauge points
elv_max = [84.2, 49.1, 54.0, 40.2, 34.9, 27.4, 21.5, 16.1, 12.9];
bot_elv = [43.9, 34.5, 30.0, 27.4, 23.1, 19.1, 11.4, 09.3, 07.5];
dep_max = elv_max - bot_elv; % the maximum depth of each gauge points

% find the maximum depth in numerical solutions
ind = 4:12;
dep_num_max = max(obj.detector.dQ(ind, :, 1), [], 2);
x = ind - 3;

markersize = 8;
linewidth = 1.5;
fontsize = 12;
figure('color', 'w'); 

% draw the maximum depth at gauge points
subplot(2,2,1); hold on;
plot(x, dep_max, 'ks-', ....
    'MarkerSize', markersize, 'LineWidth', linewidth);
plot(x, dep_num_max, 'ro--', ...
    'MarkerSize', markersize, 'LineWidth', linewidth);
legend({'Measured', 'Numerical'}, ...
    'box', 'off', 'Interpreter', 'Latex', 'FontSize', fontsize);

% draw the water depth at A, B and C
subplot(2,2,3:4); hold on;
p_h = obj.detector.draw(1, 1); set(p_h, 'Color', 'b', 'Marker', 'o');
p_h = obj.detector.draw(2, 1); set(p_h, 'Color', 'r', 'Marker', 's');
p_h = obj.detector.draw(3, 1); set(p_h, 'Color', 'g', 'Marker', '^');
legend({'A', 'B', 'C'}, ...
    'box', 'off', 'Interpreter', 'Latex', 'FontSize', fontsize);

% draw the time divergence
dep = obj.detector.dQ(1:3, :, 1); time = zeros(3,1);
ind_A = find( dep(1, :) > obj.hmin, 1 ); time(1) = obj.detector.time(ind_A);
ind_B = find( dep(1, :) > obj.hmin, 1 ); time(2) = obj.detector.time(ind_B);
ind_C = find( dep(1, :) > obj.hmin, 1 ); time(3) = obj.detector.time(ind_C);

subplot(2,2,2); hold on;
plot(1:2, diff(time), ...
    'ks', 'MarkerFaceColor', 'k', 'MarkerSize', markersize);
plot(1:2, [1140, 1320], ...
    'ks', 'MarkerFaceColor', 'k', 'MarkerSize', markersize);
end

