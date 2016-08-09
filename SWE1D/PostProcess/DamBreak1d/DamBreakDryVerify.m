function DamBreakDryVerify
filename = 'SWE1D.nc';
x = ncread(filename, 'x');
bot = zeros(size(x));

time = ncread(filename, 'time');

FinalTime = [4, 8, 12];

colorstr = {'b', 'g', 'r'};
linestr = {'', ''};
markerstr = {'o', 'x', 's'};

fig(1) = figure('Position', [519   578   534   189]);
fig(2) = figure('Position', [519   578   534   189]);
fig(3) = figure('Position', [519   578   534   189]);
% allocate line handle
p_h = zeros(numel(FinalTime)+1,1);

for itime = 1:numel(FinalTime)
    [~, index] = min(abs(time - FinalTime(itime)));
    fprintf('time deviation:%f\n', time(index) - FinalTime(itime));
    h = ncread(filename, 'h', [1,1,index],[inf,inf,1]);
    q = ncread(filename, 'q', [1,1,index],[inf,inf,1]);
    
    load(['DamBreakDry', num2str(FinalTime(itime)), 'mat']);
    
    % draw picture
    figure(fig(1));
    p_h(itime) = plot(x(:), h(:), [colorstr{itime}, linestr{1}, markerstr{itime}],...
        'MarkerSize', 5);
    hold on;
    plot(x1, h1, 'k');
    plot(x, bot, 'k')
    xlabel('x', 'Interpreter', 'Latex');
    ylabel('$\eta$', 'Interpreter', 'Latex');
    

    figure(fig(2));
    plot(x, q, [colorstr{itime}, linestr{1}, markerstr{itime}],...
        'MarkerSize', 5);
    hold on;
    plot(x1, q1, 'k')
    xlabel('x', 'Interpreter', 'Latex');
    ylabel('q', 'Interpreter', 'Latex');
    ylim([0, 35]);
    
    figure(fig(3));
    u = q./h; u(h<eps) = 0;
    plot(x, u, [colorstr{itime}, linestr{2}, markerstr{itime}],...
        'MarkerSize', 5);
    hold on;
    plot(x1, u1, 'k')
    xlabel('x', 'Interpreter', 'Latex');
    ylabel('u', 'Interpreter', 'Latex');
    ylim([0, 25]);
end% for

figure(fig(1));
p_h(end) = plot(x1, h1, 'k');

str = {'t = 4s', 't = 8s', 't = 12s', 'Exact'};
t_h = legend(p_h, str);
set(t_h, 'Box', 'off');
end% func

