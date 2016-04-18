function DamBreakDryVerify
filename = 'SWE1D.nc';
x = ncread(filename, 'x');
bedElevation = zeros(size(x));

time = ncread(filename, 'time');

FinalTime = [4, 8, 12];

for itime = 1:numel(FinalTime)
    [~, index] = min(abs(time - FinalTime(itime)));
    fprintf('time deviation： %f\n', time(index) - FinalTime(itime));
    h = ncread(filename, 'h', [1, index],[inf, 1]);
    q = ncread(filename, 'q', [1, index],[inf, 1]);
    
    load(['DamBreakDry', num2str(FinalTime(itime)), 'mat']);
    
    % draw picture
    figure; subplot(3,1,1);
    plot(x, h, '-b.'); hold on;
    plot(x1, h1, 'k')
    plot(x, bedElevation, 'k')

    subplot(3,1,2);
    plot(x, q, '-b.'); hold on;
    plot(x1, q1, 'k')
    
    subplot(3,1,3);
    u = q./h; u(h<eps) = 0;
    plot(x, u, '-b.'); hold on;
    plot(x1, u1, 'k')
    
end% for
end% func