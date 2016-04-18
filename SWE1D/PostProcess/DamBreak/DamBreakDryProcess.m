function DamBreakDryProcess
filename = 'SWE1D.nc';
x = ncread(filename, 'x');
bx = ncread(filename, 'bx');
time = ncread(filename, 'time');
startStep = 1; itime = startStep;

% get results
h = ncread(filename, 'h', [1, itime],[inf, 1]);
q = ncread(filename, 'q', [1, itime],[inf, 1]);
u = q./h; u(h<=eps) = 0;
s = ncread(filename, 'status', [1, itime], [inf, 1]);
% flat bottom
bedElevation = zeros(size(x));

% draw picture
figure; subplot(4,1,1);
p_h = plot(x, h+bedElevation, '-b.'); hold on;
plot(x, bedElevation, 'k')

subplot(4,1,2);
p_q = plot(x, q, '-b.'); hold on;

subplot(4,1,3);
p_u = plot(x, u, '-b.');

subplot(4,1,4);
p_s = plot(bx, s, 'k.', 'MarkerSize', 12);

camera_on = 0;
if camera_on
    writerObj = VideoWriter('Parabolic.avi');
    writerObj.FrameRate = 5;
    open(writerObj);
end

% time advance
for itime = 1:numel(time)
    h = ncread(filename, 'h', [1, itime],[inf, 1]);
    q = ncread(filename, 'q', [1, itime],[inf, 1]);
    u = q./h; u(h<=eps) = 0;
    s = ncread(filename, 'status', [1, itime], [inf, 1]);

    set(p_s, 'YData', s);
    set(p_h, 'YData', h+bedElevation);
    set(p_q, 'YData', q);
    set(p_u, 'YData', u);
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