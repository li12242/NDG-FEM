% close all
filename = 'SWE1D.nc';

time = ncread(filename, 'time');
% timestep = numel(time)
x = ncread(filename, 'x'); nx = numel(x);
startStep = 1;
itime = startStep;
h = ncread(filename, 'h', [1, itime],[inf, 1]);
q = ncread(filename, 'q', [1, itime],[inf, 1]);

% % ParabolicBowl
% a = 600; h0 = 10;
% bedElevation = h0.*(x.^2./a^2 - 1);
bedElevation = zeros(size(x));
flag = (x >= 8) & (x <=12);
bedElevation(flag) = 0.2 - 0.05*(x(flag) -10).^2;

figure
subplot(2,1,1); p_h = plot(x, h+bedElevation, '-b.'); hold on;
plot(x, bedElevation, 'k')
subplot(2,1,2); p_q = plot(x, q, '-r.');

for itime = startStep:4:numel(time)
    h = ncread(filename, 'h', [1, itime],[inf, 1]);
    q = ncread(filename, 'q', [1, itime],[inf, 1]);
    set(p_h, 'YData', h+bedElevation);
    set(p_q, 'YData', q);
    drawnow;
    
%     if( any((h>20)) )
%         keyboard
%     end% if

    fprintf('Processing: %f ...\n', itime/numel(time))
%     if itime > 6
%         keyboard
%     end
end% for