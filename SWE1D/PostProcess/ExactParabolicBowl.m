close all
filename = 'SWE1D.nc';

time = ncread(filename, 'time');
% timestep = numel(time)
x = ncread(filename, 'x'); nx = numel(x);
startStep = 1;
itime = startStep;
h = ncread(filename, 'h', [1, itime],[inf, 1]);
q = ncread(filename, 'q', [1, itime],[inf, 1]);

% ParabolicBowl
a = 600; h0 = 10;
bedElevation = h0.*(x.^2./a^2 - 1);
T = 269; % Parabolic Bowl period

% exact solution
hDelta = 0.0;
g = 9.8; B = 5; h0 = 10; a = 600;
w = sqrt(2*g*h0)./a;
% z = zeros(size(mesh.x));
z1 = -(4*B*w).*x./(4*g);
h1 = z1 - bedElevation;
h1(h1<hDelta) = hDelta;

z2 = (4*B*w).*x./(4*g);
h2 = z2 - bedElevation;
h2(h2<hDelta) = hDelta;

z3 = 1.22*ones(size(x));
h3 = z3 - bedElevation;
h3(h3<hDelta) = hDelta;

% plot figure
exTime = [T/2, T*3/4, T];
for n = 1:numel(exTime)
    % find result
    itime = find( abs(time - exTime(n))< 0.15 );
    
    % get result
    h = ncread(filename, 'h', [1, itime(end)],[inf, 1]);
    q = ncread(filename, 'q', [1, itime(end)],[inf, 1]);
    
    % plot
    figure
    subplot(2,1,1); 
    p_h = plot(x, h+bedElevation, '-b.'); hold on;
    plot(x, bedElevation, 'k')
    switch n
        case 1
            plot(x, bedElevation+h2, 'k--')
        case 2
%             h3 = h3 + 1.22;
            plot(x, bedElevation+h3, 'k--')
        case 3
            plot(x, bedElevation+h1, 'k--')
    end
    subplot(2,1,2); p_q = plot(x, q, '-r.');
end% for