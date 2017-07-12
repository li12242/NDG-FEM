function draw_section( obj, f_Q, time, marker)
%DRAW_SECTION Summary of this function goes here
%   Detailed explanation goes here

Nintp = 200; dx = 2/Nintp;
xd = linspace(-1+dx, 1-dx, Nintp)'; 
yd = xd;
detector = ndg_utility.detector.detector2d(obj.mesh, ...
    xd, yd, obj.ftime-1e-3, obj.ftime, obj.Nfield);
detector.collect( f_Q, time );

figure(1); hold on;
plot(xd, detector.dQ(:, 1), [marker, '-']);
end

