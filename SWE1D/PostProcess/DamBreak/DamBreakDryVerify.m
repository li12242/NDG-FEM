function DamBreakDryVerify
filename = 'SWE1D.nc';
x = ncread(filename, 'x');
bedElevation = zeros(size(x));

time = ncread(filename, 'time');

FinalTime = [4, 8, 12];

colorstr = {'b', 'g', 'r'};
linestr = {'.-', '.'};

fig(1) = figure('Position', [519   578   534   189]);
fig(2) = figure('Position', [519   578   534   189]);
fig(3) = figure('Position', [519   578   534   189]);
% allocate line handle
p_h = zeros(size(FinalTime));
for itime = 1:numel(FinalTime)
    [~, index] = min(abs(time - FinalTime(itime)));
    fprintf('time deviation： %f\n', time(index) - FinalTime(itime));
    h = ncread(filename, 'h', [1, index],[inf, 1]);
    q = ncread(filename, 'q', [1, index],[inf, 1]);
    
    load(['DamBreakDry', num2str(FinalTime(itime)), 'mat']);
    
    % draw picture
    figure(fig(1));
    p_h(itime) = plot(x, h, [colorstr{itime}, linestr{1}]); hold on;
    plot(x1, h1, 'k')
    plot(x, bedElevation, 'k')
    xlabel('x', 'Interpreter', 'Latex');
    ylabel('$\eta$', 'Interpreter', 'Latex');

    figure(fig(2));
    plot(x, q, [colorstr{itime}, linestr{1}]); hold on;
    plot(x1, q1, 'k')
    xlabel('x', 'Interpreter', 'Latex');
    ylabel('q', 'Interpreter', 'Latex');
    
    figure(fig(3));
    u = q./h; u(h<eps) = 0;
    plot(x, u, [colorstr{itime}, linestr{2}]); hold on;
    plot(x1, u1, 'k')
    xlabel('x', 'Interpreter', 'Latex');
    ylabel('u', 'Interpreter', 'Latex');
    
end% for

str = {'t = 4s', 't = 8s', 't = 12s'};
t_h = legend(p_h, str);
set(t_h, 'Box', 'off');
end% func

