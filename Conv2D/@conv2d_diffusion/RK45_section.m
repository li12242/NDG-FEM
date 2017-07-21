function [ obj ] = RK45_section( obj )
%RK45_SECTION 采用 RK45 求解时间离散并 
%   Detailed explanation goes here

rk4a = [            0.0 ...
        -567301805773.0/1357537059087.0 ...
        -2404267990393.0/2016746695238.0 ...
        -3550918686646.0/2091501179385.0  ...
        -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
         5161836677717.0/13612068292357.0 ...
         1720146321549.0/2090206949498.0  ...
         3134564353537.0/4481467310338.0  ...
         2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
         1432997174477.0/9575080441755.0 ...
         2526269341429.0/6820363962896.0 ...
         2006345519317.0/3224310063776.0 ...
         2802321613138.0/2924317926251.0];

time = 0;
ftime = obj.ftime;
f_Q = obj.f_Q;
dt = time_interval( obj, f_Q );
resQ = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);

ticktime = [0, 0.5, 1, 1.5, 1.99];
marker = 'ox+*s';
if numel(marker) < numel(ticktime) % check there are enough markers.
    error('please define more markers to show multiple time step.');
end
flag = true(size(ticktime));

while(time < ftime)
    if(time + dt > ftime)
        dt = ftime - time;
    end
    for INTRK = 1:5
        %tloc = time + rk4c(INTRK)*dt;
        rhsQ = rhs_term(obj, f_Q);
        resQ = rk4a(INTRK)*resQ + dt*rhsQ;
        
        f_Q = f_Q + rk4b(INTRK)*resQ;
    end
    time = time + dt;
    if any(flag)
        tind = find(flag, 1);
        if(time > ticktime(tind))
            draw_snapshot(obj, f_Q, ['time_', num2str(tind)]);
            draw_section(obj, f_Q, time, marker(tind));
            flag(tind) = false;
        end
    end
    %obj.draw(f_Q); drawnow; 
end

obj.f_Q = f_Q;
section_plot_legend();
end

function dt = time_interval(obj, f_Q)
spe = obj.character_len(f_Q); % Jacobian characteristic length
dt = bsxfun(@times, sqrt(obj.mesh.vol)/(2*obj.mesh.cell.N+1), 1./spe);
dt = min( min( dt ) )*0.8;
end% func

function draw_snapshot( obj, f_Q, pic_name)
fig = figure(); % new figure
obj.draw( f_Q );
grid on; box on; 
view([18, 45]); zlim([-0.2, 1]);
xlabel('$x$', 'Fontsize', 16, 'Interpreter', 'Latex');
ylabel('$y$', 'Fontsize', 16, 'Interpreter', 'Latex');
zlabel('$C$', 'Fontsize', 16, 'Interpreter', 'Latex');
print(fig, pic_name, '-r300', '-dtiff')
close(fig);
end

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

function section_plot_legend()
grid on; box on;
xlabel('$x$', 'Fontsize', 16, 'Interpreter', 'Latex');
ylabel('$C$', 'Fontsize', 16, 'Interpreter', 'Latex');
legend({'$t = 0 s$', '$t = 0.5 s$', ...
    '$t = 1 s$', '$t = 1.5 s$', '$t = 2 s$'}, ...
    'Location', 'eastoutside', ...
    'Interpreter', 'Latex',...
    'FontSize', 16);
end
