%% errors on square meshes

close all

index = 5;
data = [40	160	1.77e-02	0.00	1.85e+01	0.00
60	240	8.83e-03	1.72	1.79e+01	0.08
80	320	4.72e-03	2.18	1.58e+01	0.43
100	400	2.73e-03	2.46	1.36e+01	0.67];

n = 1;
dx = sqrt(data(:, 1).^2*(n+1).^2);
loglog(dx, data(:, index), '+-', 'Markersize', 5); hold on;

data = [40	360	1.01e-03	0.00	1.92e+00	0.00
60	540	1.96e-04	4.06	7.82e-01	2.22
80	720	6.69e-05	3.73	4.66e-01	1.80
100	900	3.16e-05	3.36	3.34e-01	1.49];

n = n+ 1;
dx = sqrt(data(:, 1).^2*(n+1).^2);
loglog(dx, data(:, index), 'o-', 'Markersize', 5);

data = [40	640	4.99e-05	0.00	1.46e-01	0.00
60	960	8.73e-06	4.30	5.30e-02	2.51
80	1280	2.68e-06	4.10	2.90e-02	2.09
100	1600	1.09e-06	4.05	1.81e-02	2.11];

n = n+ 1;
dx = sqrt(data(:, 1).^2*(n+1).^2);
loglog(dx, data(:, index), '<-', 'Markersize', 5);

data = [40	1000	3.49e-06	0.00	1.59e-02	0.00
60	1500	4.45e-07	5.08	4.43e-03	3.16
80	2000	1.05e-07	5.02	1.84e-03	3.05
100	2500	3.44e-08	4.99	9.37e-04	3.02];

n = n+ 1;
dx = sqrt(data(:, 1).^2*(n+1).^2);
loglog(dx, data(:, index), '*-', 'Markersize', 5);

data = [40	1440	7.80e-07	0.00	3.15e-03	0.00
60	2160	2.32e-08	8.67	3.28e-04	5.58
80	2880	4.08e-09	6.04	9.80e-05	4.20
100	3600	1.08e-09	5.97	4.04e-05	3.97];

n = n+ 1;
dx = sqrt(data(:, 1).^2*(n+1).^2);
loglog(dx, data(:, index), '^-', 'Markersize', 5);

data = [40	1960	8.48e-07	0.00	2.52e-03	0.00
60	2940	7.39e-07	0.34	3.40e-03	-0.74
80	3920	6.85e-08	8.26	6.05e-04	6.00
100	4900	1.19e-08	7.85	1.15e-04	7.45];

n = n+ 1;
dx = sqrt(data(:, 1).^2*(n+1).^2);
loglog(dx, data(:, index), 'd-', 'Markersize', 5);

t = legend('p = 1', 'p = 2', 'p = 3', 'p = 4', 'p = 5', 'p = 6');
set(t, 'box', 'off')

xlabel('$\sqrt{DOFs}$', 'Interpreter', 'Latex')
ylabel('$L_{\infty}$', 'Interpreter', 'Latex')

%% errors on triangle

index = 5;

data = [40 	 9600	 1.56e-02 	0.00 	 2.43e+01 	0.00 	
60 	 21600  7.59e-03 	1.78 	 2.36e+01 	0.07 	
80 	 38400  4.08e-03  2.16 	 2.10e+01 	0.39 	
100  60000  2.41e-03 	2.36 	 1.87e+01 	0.52];

dx = sqrt(data(:, 2));
loglog(dx, data(:, index), 'b+-', 'Markersize', 5); hold on;

data = [40 	 19200  8.27e-02 	0.00 	 1.97e+02 	0.00 
60 	 43200  8.27e-02 	-0.00 	 4.43e+02 	-2.00 	
80 	 76800  9.98e+14 	-128.71 5.51e+17 	-120.82
100  120000 4.35e+19 	-47.87  3.02e+22 	-48.90 ];

% compared square meshes
data = [40	160	1.77e-02	0.00	1.85e+01	0.00
60	240	8.83e-03	1.72	1.79e+01	0.08
80	320	4.72e-03	2.18	1.58e+01	0.43
100	400	2.73e-03	2.46	1.36e+01	0.67];

n = 1;
dx = sqrt(data(:, 1).^2*(n+1).^2);
loglog(dx, data(:, index), 'r+-', 'Markersize', 5); 

data = [40	360	1.01e-03	0.00	1.92e+00	0.00
60	540	1.96e-04	4.06	7.82e-01	2.22
80	720	6.69e-05	3.73	4.66e-01	1.80
100	900	3.16e-05	3.36	3.34e-01	1.49];

n = n+ 1;
dx = sqrt(data(:, 1).^2*(n+1).^2);
loglog(dx, data(:, index), 'ro-', 'Markersize', 5);

t = legend('triangle p = 1', 'square p = 1', 'square p = 2');
set(t, 'box', 'off')
xlabel('$\sqrt{DOFs}$', 'Interpreter', 'Latex')
% ylabel('$L_2$', 'Interpreter', 'Latex')
ylabel('$L_{\infty}$', 'Interpreter', 'Latex')