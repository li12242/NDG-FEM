function TesunamiRunupVerify
SingleResult;
end

function SingleResult
filename = 'SWE1DTsunamiRunup.nc';
x = ncread(filename, 'x'); x = reshape(x, 2, numel(x)/2);
bedElevation = -0.1*x;
time = ncread(filename, 'time');

FinalTime = [160, 175, 220];

for itime = 1:numel(FinalTime)
    [~, index] = min(abs(time - FinalTime(itime)));
    fprintf('time deviation %f\n', time(index) - FinalTime(itime));
    h = ncread(filename, 'h', [1, index],[inf, 1]);
    q = ncread(filename, 'q', [1, index],[inf, 1]);
    h = reshape(h, 2, numel(x)/2); q = reshape(q, 2, numel(x)/2);
    eta = h + bedElevation; u = q./h; u(h<1e-2) = 0; u(u>50) = 0;
    
    load(['t', num2str(FinalTime(itime)), '.mat']);
    
    % draw picture
    figure('Position', [266   506   556   252]);
    bed1 = -0.1*x1; eta1(eta1 < bed1) = bed1(eta1 < bed1);
    plot(x1, eta1, 'k.', 'MarkerSize', 8); hold on;
    plot(x, eta, '-r.', 'LineWidth', 1);
    plot(x, bedElevation, 'k')
    xlim([-400, 800]);
    t = legend('Exact', 'HREF');
    set(t, 'box', 'off');
        
    xlabel('x (m)', 'Interpreter', 'Latex');
    ylabel('$\eta$ (m)', 'Interpreter', 'Latex');
    
    figure('Position', [266   506   556   252]);
    plot(x1, u1, 'k.', 'MarkerSize', 8); hold on;
    plot(x, u, '-r.', 'LineWidth', 1); 
    ylim([-20, 15])
    xlim([-400, 800]);
    xlabel('x (m)', 'Interpreter', 'Latex');
    ylabel('u (m/s)', 'Interpreter', 'Latex');
    
end% for
end

function ComparedResult
filename = 'SWE1D.nc';
filename1 = 'SWE1DTsunamiRunup1.nc';
x = ncread(filename, 'x');
bedElevation = -0.1*x;

time = ncread(filename, 'time');
time1 = ncread(filename1, 'time');

FinalTime = [160, 175, 220];

for itime = 1:numel(FinalTime)
    [~, index] = min(abs(time - FinalTime(itime)));
    fprintf('time deviation %f\n', time(index) - FinalTime(itime));
    h = ncread(filename, 'h', [1, index],[inf, 1]);
    q = ncread(filename, 'q', [1, index],[inf, 1]);
    eta = h + bedElevation; u = q./h; u(h<1e-3) = 0; u(u>50) = 0;
    
    [~, index1] = min(abs(time1 - FinalTime(itime)));
    fprintf('time deviation %f\n', time1(index1) - FinalTime(itime));
    hr = ncread(filename1, 'h', [1, index1],[inf, 1]);
    qr = ncread(filename1, 'q', [1, index1],[inf, 1]);
    etar = h + bedElevation; ur = qr./hr; ur(hr<1e-3) = 0; ur(ur>50) = 0;
    
    load(['t', num2str(FinalTime(itime)), '.mat']);
    
    % draw picture
    figure; subplot(2,1,1);
    
    plot(x, eta, '-bx'); hold on;
    plot(x, etar, '-ro');
    plot(x1, eta1, 'k'); 
    plot(x, bedElevation, 'k')
    xlim([-400, 800]);
    legend('convectional wet/dry treatment', 'h-refined wet/dry treatment', ...
        'analytical solution')
    xlabel('x (m)', 'Interpreter', 'Latex');
    ylabel('$\eta$ (m)', 'Interpreter', 'Latex');
    
    subplot(2,1,2);
    plot(x1, u1, 'k'); hold on;
    plot(x, u, '-bx'); 
    plot(x, ur, '-ro')
    
    xlim([-400, 800]);
    xlabel('x (m)', 'Interpreter', 'Latex');
    ylabel('u (m/s)', 'Interpreter', 'Latex');
    
end% for
end% func