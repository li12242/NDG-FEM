function LakeAtRestPostprocess
filename = 'SWE1D_LakeAtRest.nc';
filename1 = 'SWE1D_LakeAtRest1.nc';
time = ncread(filename, 'time');
% timestep = numel(time)
x = ncread(filename, 'x');
startStep = 1;
itime = startStep;

h = ncread(filename, 'h', [1, itime],[inf, 1]);
q = ncread(filename, 'q', [1, itime],[inf, 1]);
h1 = ncread(filename1, 'h', [1, itime],[inf, 1]);
q1 = ncread(filename1, 'q', [1, itime],[inf, 1]);

% bottom topography
bedElevation = zeros(size(x));
a = 1.2; rm = 0.4; r = abs(x - 0.5);
index = (r < rm);
bedElevation(index) = a*exp(-0.5./(rm.^2 - r(index).^2))./exp(-0.5./rm^2);

figure; subplot(2,1,1);
p_h = plot(x, h+bedElevation, '-b.'); hold on;
p_h1 = plot(x, h1+bedElevation, '-r.');
plot(x, bedElevation, 'k')

subplot(2,1,2);
p_q = plot(x, q, '-b.'); hold on;
p_q1 = plot(x, q1, '-r.'); hold on;
camera_on = 0;
if camera_on
    writerObj = VideoWriter('Parabolic.avi');
    writerObj.FrameRate = 5;
    open(writerObj);
end

for itime = 1:numel(time)
    h = ncread(filename, 'h', [1, itime],[inf, 1]);
    q = ncread(filename, 'q', [1, itime],[inf, 1]);
    h1 = ncread(filename1, 'h', [1, itime],[inf, 1]);
    q1 = ncread(filename1, 'q', [1, itime],[inf, 1]);
    
    set(p_h, 'YData', h+bedElevation);
    set(p_q, 'YData', q);
    set(p_h1, 'YData', h1+bedElevation);
    set(p_q1, 'YData', q1./100);
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
end% func