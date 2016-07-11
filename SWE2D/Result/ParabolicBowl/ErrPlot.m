% Triangle
triData_H = [80 	 25600 	 5.36e-03 	 0.00 	 5.68e-02 	 0.00
100 	 40000 	 4.02e-03 	 1.29 	 4.78e-02 	 0.78
120 	 57600 	 3.21e-03 	 1.23 	 4.04e-02 	 0.92
160 	 102400 	 2.40e-03 	 1.01 	 2.99e-02 	 1.04
200 	 160000 	 2.06e-03 	 0.69 	 2.37e-02 	 1.05];

triData_Q = [80 	 25600 	 1.06e-02 	 0.00 	 6.80e-02 	 0.00
100 	 40000 	 8.69e-03 	 0.88 	 5.85e-02 	 0.67
120 	 57600 	 7.48e-03 	 0.82 	 5.30e-02 	 0.55
160 	 102400 	 6.31e-03 	 0.59 	 4.66e-02 	 0.45
200 	 160000 	 5.55e-03 	 0.57 	 4.14e-02 	 0.52];

% Quadrilaterial
quadData_H = [80 	 25600 	 4.54e-03 	 0.00 	 4.01e-02 	 0.00
100 	 40000 	 3.36e-03 	 1.35 	 3.62e-02 	 0.46
120 	 57600 	 2.53e-03 	 1.55 	 2.50e-02 	 2.03
160 	 102400 	 1.84e-03 	 1.11 	 2.23e-02 	 0.39
200 	 160000 	 1.50e-03 	 0.92 	 1.85e-02 	 0.85];

quadData_Q = [80 	 25600 	 1.01e-02 	 0.00 	 4.61e-02 	 0.00
100 	 40000 	 7.45e-03 	 1.36 	 3.50e-02 	 1.23
120 	 57600 	 5.83e-03 	 1.35 	 2.77e-02 	 1.29
160 	 102400 	 4.30e-03 	 1.06 	 2.16e-02 	 0.86
200 	 160000 	 3.52e-03 	 0.90 	 1.77e-02 	 0.90];

index = [3,5];
n     = 1;
yStr  = {'$L_2$', '$L_{\infty}$'};
tStr  = {'$\rm{H}$', '$\rm{Q_x}$'};
% Plot H
for i = 1:2
    figure
    dx = sqrt(triData_H(:, 1).^2*2 * (n+1)*(n+2)/2);
    loglog(dx, triData_H(:, index(i)), 'bo-',...
        'Markersize', 8, 'MarkerFaceColor', 'b'); hold on;
    
    loglog(dx, triData_Q(:, index(i)), 'b+-',...
        'Markersize', 8, 'MarkerFaceColor', 'b');

    dx = sqrt(quadData_H(:, 1).^2*(n+1).^2);
    loglog(dx, quadData_H(:, index(i)), 'ro-', ...
        'Markersize', 8, 'MarkerFaceColor', 'r');
    
    loglog(dx, quadData_Q(:, index(i)), 'r+-', ...
        'Markersize', 8, 'MarkerFaceColor', 'r');
    
    t = legend('Tri $h$', 'Tri $q_x$', 'Quad $h$', 'Quad $q_x$', ...
        'Location', 'NorthEast');
    set(t, 'box', 'off', 'Interpreter', 'Latex')
    xlabel('$\sqrt{DOFs}$', 'Interpreter', 'Latex')
    ylabel(yStr{i}, 'Interpreter', 'Latex')
    grid on;
    title(tStr{i},  'Interpreter', 'Latex')
end