function result( obj )
%RESULT Summary of this function goes here
%   Detailed explanation goes here

len = 1000;
Nintp = 300; 
dx = len/(Nintp+1);
xd = linspace(dx, len-dx, Nintp)'; 
yd = ones(Nintp, 1)*10;
T = [cos(obj.theta), sin(obj.theta); -sin(obj.theta), cos(obj.theta)];
temp = T*[ xd(:)'; yd(:)' ];
xd(:) = temp(1, :);
yd(:) = temp(2, :);

detector = ndg_utility.detector.detector2d(obj.mesh, ...
    xd, yd, obj.ftime/2, obj.ftime, obj.Nfield);
detector.collect( obj.f_Q, obj.ftime );

% get the exact solutions
[ f_ext ] = obj.ext_func( obj.ftime );
detector.collect( f_ext, obj.ftime );

h_num = detector.dQ(:, 1, 1); hu_num = detector.dQ(:, 1, 2);
h_ext = detector.dQ(:, 2, 1); hu_ext = detector.dQ(:, 2, 2);

figure('color', 'w'); 
subplot(2,1,1); hold on;
plot(xd, h_ext, 'r-', 'LineWidth', 2);
plot(xd, h_num, 'bo-', 'MarkerSize', 6);
box on; grid on;

subplot(2,1,2); hold on;
plot(xd, hu_ext, 'r-', 'LineWidth', 2);
plot(xd, hu_num, 'bo-', 'MarkerSize', 6);
box on; grid on;
end
