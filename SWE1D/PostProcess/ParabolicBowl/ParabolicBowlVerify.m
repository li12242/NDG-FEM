function ParabolicBowlVerify

filename = 'SWE1D_2_100.nc';
time = ncread(filename, 'time');
x = ncread(filename, 'x');

% bottom topography
g = 9.8; B = 5; h0 = 10; a = 3000; w = sqrt(2*g*h0)./a;
bedElevation = h0.*(x.^2./a^2 - 1);

T = 269;

FinalTime = [T/4, T/2, T];

for itime = 1:numel(FinalTime)
    [~, index] = min(abs(time - FinalTime(itime)));
    
    h = ncread(filename, 'h', [1, index],[inf, 1]);
    q = ncread(filename, 'q', [1, index],[inf, 1]);
    u = q./h; u(h<1e-3) = 0.0;
    
    % get exact solution
    [he, qe, ue] = exactParabolicBowlSolution(FinalTime(itime), x, bedElevation);
    
    % draw picture
    figure; subplot(3,1,1)
    plot(x, he+bedElevation, 'r.-'); hold on;
    plot(x, bedElevation, 'k')
    plot(x, h+bedElevation, 'b'); 
    xlabel('x', 'Interpreter', 'Latex');
    ylabel('$\eta$', 'Interpreter', 'Latex');

    subplot(3,1,2)
    plot(x, q, 'b'); hold on;
    plot(x, qe, 'r');
    xlabel('x', 'Interpreter', 'Latex');
    ylabel('q', 'Interpreter', 'Latex');
    
    subplot(3,1,3)
    plot(x, u, 'b'); hold on;
    plot(x, ue, 'r');
    xlabel('x', 'Interpreter', 'Latex');
    ylabel('u', 'Interpreter', 'Latex');
end% for
end% func
