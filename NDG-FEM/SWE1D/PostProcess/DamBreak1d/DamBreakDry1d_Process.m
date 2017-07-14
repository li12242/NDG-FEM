function DamBreakDry1d_Process
filename = 'SWE1D_DamBreakDry_1000.nc';
x = ncread(filename, 'x');
bot = ncread(filename, 'bot');
time = ncread(filename, 'time');
startStep = 1; itime = startStep;

% get results
h = ncread(filename, 'h', [1, 1, itime],[inf, inf, 1]);
q = ncread(filename, 'q', [1, 1, itime],[inf, inf, 1]);
u = q./h; u(h<=eps) = 0;

% draw picture
figure; 
subplot(3,1,1);
p_h = plot(x(:), h(:)+bot(:), '-b.'); hold on;
plot(x(:), bot(:), 'k')
% ylim([99.8, 101])

subplot(3,1,2);
p_q = plot(x(:), q(:), '-b.'); hold on;
% xlim([0, 1e3]); ylim([0, 30]);

subplot(3,1,3);
p_u = plot(x(:), u(:), '-b.');
% xlim([0, 1e3]); ylim([0, 20]);

camera_on = 0;
if camera_on
    writerObj = VideoWriter('Parabolic.avi');
    writerObj.FrameRate = 5;
    open(writerObj);
end

% time advance
for itime = 1:numel(time)
    h = ncread(filename, 'h', [1,1,itime],[inf,inf,1]);
    if any(any(isnan(h)))
        keyboard;
    end% if
    q = ncread(filename, 'q', [1,1,itime],[inf,inf,1]);
    u = q./h; u(h<=eps) = 0;
    
    set(p_h, 'YData', h(:)+bot(:));
    set(p_q, 'YData', q(:));
    set(p_u, 'YData', u(:));
    drawnow;

    if camera_on
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
    end

    fprintf('Processing: %f ...\n', itime/numel(time))
end% for

if camera_on
    close(writerObj);
end
end% func