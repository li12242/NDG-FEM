function result( obj )
%RESULT Summary of this function goes here
%   Detailed explanation goes here

% 读取实测数据
file = 'SWE2D/@swe2d_tsuami/mesh/output_ch5-7-9.xls';
data = xlsread(file);
data(:, [2,3,4]) = data(:, [2,3,4])./100; % 转换为 m

markersize = 6;
linewidth = 3;
fontsize = 14;
% 绘制实测点 ch5
figure('color', 'w', 'Position', [100, 170, 625, 635]);
subplot(3,1,1);
plot(data(:, 1), data(:, 2), 'ro', 'MarkerSize', markersize);
hold on; grid on;
p_h = obj.detector.draw(1, 1);
set(p_h, 'LineWidth', linewidth);
legend({'Measured', 'Numerical'}, ...
    'box', 'off', 'Location', 'NorthWest', ...
    'FontSize', fontsize);
xlim([0, obj.ftime]);
% 绘制实测点 ch7
subplot(3,1,2);
plot(data(:, 1), data(:, 3), 'ro', 'MarkerSize', markersize);
hold on; grid on;
p_h = obj.detector.draw(2, 1);
set(p_h, 'LineWidth', linewidth);
xlim([0, obj.ftime]);

% 绘制实测点 ch9
subplot(3,1,3);
plot(data(:, 1), data(:, 4), 'ro', 'MarkerSize', markersize);
hold on; grid on;
p_h = obj.detector.draw(3, 1);
set(p_h, 'LineWidth', linewidth);
xlim([0, obj.ftime]);

end

