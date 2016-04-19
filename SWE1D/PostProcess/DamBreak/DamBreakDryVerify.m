function DamBreakDryVerify
filename = 'SWE1D_Dry.nc';
x = ncread(filename, 'x');
bedElevation = zeros(size(x));

time = ncread(filename, 'time');

FinalTime = [4, 8, 12];

colorstr = {'-b.', '-g.', '-r.'};
for itime = 1:numel(FinalTime)
    [~, index] = min(abs(time - FinalTime(itime)));
    fprintf('time deviation： %f\n', time(index) - FinalTime(itime));
    h = ncread(filename, 'h', [1, index],[inf, 1]);
    q = ncread(filename, 'q', [1, index],[inf, 1]);
    
    load(['DamBreakDry', num2str(FinalTime(itime)), 'mat']);
    
    % draw picture
    subplot(3,1,1);
    plot(x, h, colorstr{itime}); hold on;
    plot(x1, h1, 'k')
    plot(x, bedElevation, 'k')
    xlabel('x', 'Interpreter', 'Latex');
    ylabel('$\eta$', 'Interpreter', 'Latex');

    subplot(3,1,2);
    plot(x, q, colorstr{itime}); hold on;
    plot(x1, q1, 'k')
    xlabel('x', 'Interpreter', 'Latex');
    ylabel('q', 'Interpreter', 'Latex');
    
    subplot(3,1,3);
    u = q./h; u(h<eps) = 0;
    plot(x, u, colorstr{itime}); hold on;
    plot(x1, u1, 'k')
    xlabel('x', 'Interpreter', 'Latex');
    ylabel('u', 'Interpreter', 'Latex');
    
end% for
end% func