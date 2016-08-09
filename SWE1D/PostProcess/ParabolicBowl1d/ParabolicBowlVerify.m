function ParabolicBowlVerify

filename = 'SWE1D_ParabolicBowl.nc';
time = ncread(filename, 'time');
x = ncread(filename, 'x');

np = 50;
xe = linspace(min(x), max(x), np);

% bottom topography
g = 9.81; B = 5; h0 = 10; a = 3000; w = sqrt(2*g*h0)./a;
b = h0.*(x.^2./a^2 - 1);
be = h0.*(xe.^2./a^2 - 1);

T = 2*pi*a/sqrt(2*g*h0);

FinalTime = [T/4, T/2, T];

for itime = 1:numel(FinalTime)
    [~, index] = min(abs(time - FinalTime(itime)));
    
    h = ncread(filename, 'h', [1, index],[inf, 1]);
    q = ncread(filename, 'q', [1, index],[inf, 1]);
    u = q./h; u(h<1e-3) = 0.0;
    
    % get exact solution
    [he, qe, ue] = exactParabolicBowlSolution(FinalTime(itime), xe, be);
    
    % draw picture
    figure('Position', [771   467   554   293]);
    plot(xe, he+be, 'b-', 'LineWidth', 2); hold on;
    plot(x, h+b, 'r--', 'LineWidth', 2); 
    t = legend('Exact', 'HREF');
    plot(x, b, 'k')
    xlim([-4500, 4500])
    set(t, 'box', 'off');
    xlabel('x', 'Interpreter', 'Latex');
    ylabel('$\eta$', 'Interpreter', 'Latex');

    figure('Position', [771   467   554   293]);
    plot(xe, qe, 'b-', 'LineWidth', 2); hold on;
    plot(x, q, 'r--', 'LineWidth', 2); 
    xlim([-4500, 4500])
    ylim([-55, 55]);
    set(gca, 'YTick', -50:25:50)
    xlabel('x', 'Interpreter', 'Latex');
    ylabel('q', 'Interpreter', 'Latex');
    
    figure('Position', [771   467   554   293]);
    plot(xe, ue, 'b-', 'LineWidth', 2); hold on;
    plot(x, u, 'r--', 'LineWidth', 2); 
    xlim([-4500, 4500]); ylim([-10, 10]);
    xlabel('x', 'Interpreter', 'Latex');
    ylabel('u', 'Interpreter', 'Latex');
end% for
end% func
