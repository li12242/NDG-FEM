% flow condition, 
% 1 - subcritical flow
% 2 - supercritical flow
% 3 - transcritical flow
condition = 3;
[ex, exh, exz] = ExactFlowDump(condition);

% get result
filename = 'TranscriticalFlow.nc';

time = ncread(filename, 'time');
timestep = numel(time);
x = ncread(filename, 'x'); nx = numel(x);
itime = timestep;
h = ncread(filename, 'h', [1, itime],[inf, 1]);
q = ncread(filename, 'q', [1, itime],[inf, 1]);

% bottom elevation
bedElevation = zeros(size(x));
flag = (x >= 8) & (x <=12);
bedElevation(flag) = 0.2 - 0.05*(x(flag) -10).^2;

%% plot figure
% water level
figure
subplot(2,1,1)
plot(x(:), bedElevation(:)+h(:), 'b.-')
hold on
plot(ex(:), exz(:)+exh(:), 'k');
legend('DGM', 'Exact')
plot(x(:), bedElevation(:), 'k');
title('water depth')
% flux
subplot(2,1,2)
plot(x(:), q(:), 'r.-'); hold on;
switch condition
    case 1
%         q0 = 0.18; h0 = 0.5; 
        plot(ex, 0.18*ones(size(ex)), 'k-');
        ylim([0.16, 0.2]);
    case 2
%         q0 = 25.0567; h0 = 2.0; 
        plot(ex, 20.0567*ones(size(ex)), 'k-');
        ylim([18, 22]);
    case 3
%         q0 = 0.18; h0 = 0.33;
        plot(ex, 0.18*ones(size(ex)), 'k-');
        ylim([0, 2]);
end
legend('DGM', 'Exact')
title('flux')

