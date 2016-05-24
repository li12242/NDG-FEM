% %% errors on square meshes
% 
% close all
% 
% index = 5;
% data = [40	160	1.77e-02	0.00	1.85e+01	0.00
% 60	240	8.83e-03	1.72	1.79e+01	0.08
% 80	320	4.72e-03	2.18	1.58e+01	0.43
% 100	400	2.73e-03	2.46	1.36e+01	0.67];
% 
% n = 1;
% dx = sqrt(data(:, 1).^2*(n+1).^2);
% loglog(dx, data(:, index), '+-', 'Markersize', 5); hold on;
% 
% data = [40	360	1.01e-03	0.00	1.92e+00	0.00
% 60	540	1.96e-04	4.06	7.82e-01	2.22
% 80	720	6.69e-05	3.73	4.66e-01	1.80
% 100	900	3.16e-05	3.36	3.34e-01	1.49];
% 
% n = n + 1;
% dx = sqrt(data(:, 1).^2*(n+1).^2);
% loglog(dx, data(:, index), 'o-', 'Markersize', 5);
% 
% data = [40	640	4.99e-05	0.00	1.46e-01	0.00
% 60	960	8.73e-06	4.30	5.30e-02	2.51
% 80	1280	2.68e-06	4.10	2.90e-02	2.09
% 100	1600	1.09e-06	4.05	1.81e-02	2.11];
% 
% n = n+ 1;
% dx = sqrt(data(:, 1).^2*(n+1).^2);
% loglog(dx, data(:, index), '<-', 'Markersize', 5);
% 
% data = [40	1000	3.49e-06	0.00	1.59e-02	0.00
% 60	1500	4.45e-07	5.08	4.43e-03	3.16
% 80	2000	1.05e-07	5.02	1.84e-03	3.05
% 100	2500	3.44e-08	4.99	9.37e-04	3.02];
% 
% n = n+ 1;
% dx = sqrt(data(:, 1).^2*(n+1).^2);
% loglog(dx, data(:, index), '*-', 'Markersize', 5);
% 
% data = [40	1440	7.80e-07	0.00	3.15e-03	0.00
% 60	2160	2.32e-08	8.67	3.28e-04	5.58
% 80	2880	4.08e-09	6.04	9.80e-05	4.20
% 100	3600	1.08e-09	5.97	4.04e-05	3.97];
% 
% n = n+ 1;
% dx = sqrt(data(:, 1).^2*(n+1).^2);
% loglog(dx, data(:, index), '^-', 'Markersize', 5);
% 
% data = [40	1960	8.48e-07	0.00	2.52e-03	0.00
% 60	2940	7.39e-07	0.34	3.40e-03	-0.74
% 80	3920	6.85e-08	8.26	6.05e-04	6.00
% 100	4900	1.19e-08	7.85	1.15e-04	7.45];
% 
% n = n+ 1;
% dx = sqrt(data(:, 1).^2*(n+1).^2);
% loglog(dx, data(:, index), 'd-', 'Markersize', 5);
% 
% t = legend('p = 1', 'p = 2', 'p = 3', 'p = 4', 'p = 5', 'p = 6');
% set(t, 'box', 'off')
% 
% xlabel('$\sqrt{DOFs}$', 'Interpreter', 'Latex')
% ylabel('$L_{\infty}$', 'Interpreter', 'Latex')

%% errors on triangle
index = 5;
figure

data = [20 	 2400 	 4.03e-02 	 0.00 	 7.60e-01 	 0.00
40 	 9600 	 2.48e-02 	 0.70 	 5.35e-01 	 0.51
60 	 21600 	 2.05e-02 	 0.46 	 4.20e-01 	 0.60
80 	 38400 	 1.62e-02 	 0.82 	 3.28e-01 	 0.86
100 	 60000 	 1.40e-02 	 0.67 	 2.66e-01 	 0.95];

n = 1;
dx = sqrt(data(:, 1).^2*2 * (n+1)*(n+2)/2);
loglog(dx, data(:, index), 'b+-', 'Markersize', 8, 'MarkerFaceColor', 'b'); hold on;

data = [20 	 4800 	 1.41e-02 	 0.00 	 3.32e-01 	 0.00
40 	 19200 	 5.87e-03 	 0.84 	 1.82e-01 	 0.87
60 	 43200 	 2.18e-03 	 1.56 	 9.98e-02 	 1.48
80 	 76800 	 0.93e-03 	 10.68 	 6.39e-02 	 6.85
100 	 120000 	 5.18e-04 	 -4.42 	 4.24e-02 	 -4.99];

n = 2;
dx = sqrt(data(:, 1).^2*2 * (n+1)*(n+2)/2);
loglog(dx, data(:, index), 'bo-', 'Markersize', 5, 'MarkerFaceColor', 'b');

data = [20 	 8000 	 6.33e-03 	 0.00 	 3.32e-01 	 0.00
40 	 32000 	 1.19e-04 	 5.73 	 4.29e-03 	 6.27
60 	 72000 	 2.32e-05 	 4.04 	 9.25e-04 	 3.78
80 	 128000 	 7.57e-06 	 3.88 	 2.99e-04 	 3.92
100 	 200000 	 3.15e-06 	 -1.85 	 1.54e-04 	 -7.35];

n = 3;
dx = sqrt(data(:, 1).^2*2 * (n+1)*(n+2)/2);
loglog(dx, data(:, index), 'b*-', 'Markersize', 8, 'MarkerFaceColor', 'b');

% compared square meshes
data = [20 	 1600 	 3.83e-02 	 0.00 	 7.17e-01 	 0.00
40 	 6400 	 1.65e-02 	 1.21 	 3.96e-01 	 0.86
60 	 14400 	 8.30e-03 	 1.69 	 2.20e-01 	 1.45
80 	 25600 	 4.65e-03 	 2.02 	 1.28e-01 	 1.87
100 	 40000 	 2.73e-03 	 2.39 	 7.86e-02 	 2.20];

n = 1;
dx = sqrt(data(:, 1).^2*(n+1).^2);
loglog(dx, data(:, index), 'r+-', 'Markersize', 8, 'MarkerFaceColor', 'r');

data = [20 	 3600 	 9.43e-03 	 0.00 	 2.21e-01 	 0.00
40 	 14400 	 1.01e-03 	 3.22 	 2.65e-02 	 3.06
60 	 32400 	 1.96e-04 	 4.06 	 4.78e-03 	 4.23
80 	 57600 	 6.69e-05 	 3.73 	 1.68e-03 	 3.64
100 	 90000 	 3.16e-05 	 3.36 	 9.11e-04 	 2.73];

n = n+ 1;
dx = sqrt(data(:, 1).^2*(n+1).^2);
loglog(dx, data(:, index), 'ro-', 'Markersize', 5, 'MarkerFaceColor', 'r');

data = [20 	 6400 	 1.66e-03 	 0.00 	 3.21e-02 	 0.00
40 	 25600 	 4.99e-05 	 5.06 	 1.97e-03 	 4.03
60 	 57600 	 8.73e-06 	 4.30 	 4.20e-04 	 3.81
80 	 102400 	 2.68e-06 	 4.10 	 1.37e-04 	 3.88
100 	 160000 	 1.09e-06 	 4.05 	 5.73e-05 	 3.92];

n = n+ 1;
dx = sqrt(data(:, 1).^2*(n+1).^2);
loglog(dx, data(:, index), 'r*-', 'Markersize', 8, 'MarkerFaceColor', 'r');

t = legend('triangle p = 1', 'triangle p = 2', 'triangle p = 3',...
    'square p = 1', 'square p = 2', 'square p = 3', 'Location', 'SouthWest');
set(t, 'box', 'off')
xlabel('$\sqrt{DOFs}$', 'Interpreter', 'Latex')
% ylabel('$L_2$', 'Interpreter', 'Latex')
ylabel('$L_{\infty}$', 'Interpreter', 'Latex')
grid on;