function ParabolicBowlPostProcess

filename = 'SWE1D.nc';
time = ncread(filename, 'time'); x = ncread(filename, 'x');
startStep = 1; itime = startStep;

% bottom topography
a = 3000; h0 = 10; g = 9.81; B = 5; w = sqrt(2*g*h0)./a;
b = h0.*(x.^2./a^2 - 1);

% get result
h = ncread(filename, 'h', [1, itime],[inf, 1]);
q = ncread(filename, 'q', [1, itime],[inf, 1]);

np = 100;
xe = linspace(min(x), max(x), np);
be = h0.*(xe.^2./a^2 - 1);
% get exact solution
[he, qe, ue] = exactParabolicBowlSolution(time(1), xe, be);

% draw
figure; subplot(3,1,1);
p_h = plot(x, h+b, '-b'); hold on;
p_h1 = plot(xe, he+be, 'r.');
plot(x, b, 'k')

subplot(3,1,2);
p_q = plot(x, q, '-b'); hold on;
p_q1 = plot(xe, qe, 'r.');

u = q./h; u(h<=0) = 0;
subplot(3,1,3);
p_u = plot(x, u, '-b'); hold on;
p_u1 = plot(xe, ue, 'r.');
ymax = B*a*w./sqrt(2*h0*g);
ylim([-2*ymax, 2*ymax]);

camera_on = 0;
if camera_on
    writerObj = VideoWriter('Parabolic.avi');
    writerObj.FrameRate = 5;
    open(writerObj);
end

for itime = 1:numel(time)
    h = ncread(filename, 'h', [1, itime],[inf, 1]);
    q = ncread(filename, 'q', [1, itime],[inf, 1]);
    [he, qe, ue] = exactParabolicBowlSolution(time(itime), xe, be);
    
    u = q./h; u(h<=1e-3) = 0;
    set(p_h, 'YData', h+b);
    set(p_q, 'YData', q);
    set(p_u, 'YData', u);
    set(p_h1, 'YData', he+be);
    set(p_q1, 'YData', qe);
    set(p_u1, 'YData', ue);
    drawnow;
    
%     if( any((h>20)) )
%         keyboard
%     end% if

    if camera_on
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
    end

    fprintf('Processing: %f ...\n', itime/numel(time))
    
%     if itime/numel(time)> .969188
%         keyboard
%     end
%     if itime > 695
%         keyboard
%     end
end% for

if camera_on
    close(writerObj);
end

end% func