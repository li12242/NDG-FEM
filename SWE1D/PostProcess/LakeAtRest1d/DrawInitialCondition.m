function DrawInitialCondition
filename = 'SWE1D_LakeAtRest.nc';
% timestep = numel(time)
x = ncread(filename, 'x');
startStep = 1;
itime = startStep;

h = ncread(filename, 'h', [1, itime],[inf, 1]);

% bottom topography
bedElevation = zeros(size(x));
a = 1.2; rm = 0.4; r = abs(x - 0.5);
index = (r < rm);
bedElevation(index) = a*exp(-0.5./(rm.^2 - r(index).^2))./exp(-0.5./rm^2);

figure('Position',[327   558   574   187]);
p_h = plot(x, h+bedElevation, '-b.'); hold on;
plot(x, bedElevation, 'k')
ylim([-.5, 1.5])
xlabel('x'); ylabel('$\eta$');

end% func