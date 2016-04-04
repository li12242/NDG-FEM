% close all
filename = 'SWE1D.nc';
filename2 = 'SWE1D_Hrefined.nc';

time = ncread(filename2, 'time');
% timestep = numel(time)
x = ncread(filename, 'x'); nx = numel(x);
startStep = 1;
itime = startStep;
h = ncread(filename, 'h', [1, itime],[inf, 1]);
q = ncread(filename, 'q', [1, itime],[inf, 1]);

h1 = ncread(filename2, 'h', [1, itime],[inf, 1]);
q1 = ncread(filename2, 'q', [1, itime],[inf, 1]);

% % ParabolicBowl
a = 600; h0 = 10;
bedElevation = h0.*(x.^2./a^2 - 1);
% bedElevation = zeros(size(x));
% flag = (x >= 8) & (x <=12);
% bedElevation(flag) = 0.2 - 0.05*(x(flag) -10).^2;

figure
subplot(2,1,1); p_h = plot(x, h+bedElevation, '-b.'); hold on;
p_h2 = plot(x, h1+bedElevation, 'r.-');
legend('no Wet/Dry reconstruction', 'h-refined Wet/Dry reconstruction')
plot(x, bedElevation, 'k')

subplot(2,1,2); 
p_q = plot(x, q, '-b.'); hold on;
p_q2 = plot(x, q1, 'r.-');
legend('no Wet/Dry reconstruction', 'h-refined Wet/Dry reconstruction')
ylim([-1, 1])

camera_on = 1;
if camera_on
    writerObj = VideoWriter('Parabolic.avi');
    writerObj.FrameRate = 5;
    open(writerObj);
end

for itime = 1:numel(time)
    h = ncread(filename, 'h', [1, itime],[inf, 1]);
    q = ncread(filename, 'q', [1, itime],[inf, 1]);
    
    h1 = ncread(filename2, 'h', [1, itime],[inf, 1]);
    q1 = ncread(filename2, 'q', [1, itime],[inf, 1]);

    set(p_h, 'YData', h+bedElevation);
    set(p_q, 'YData', q);
    
    set(p_h2, 'YData', h1+bedElevation);
    set(p_q2, 'YData', q1);
    drawnow;
    
%     if( any((h>20)) )
%         keyboard
%     end% if

    if camera_on
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
    end

    fprintf('Processing: %f ...\n', itime/numel(time))
%     if itime > 6
%         keyboard
%     end
end% for

if camera_on
    close(writerObj);
end
% [x, eta] = ExactFlowDump(3);
% subplot(2,1,1);
% plot(x, eta, 'k');