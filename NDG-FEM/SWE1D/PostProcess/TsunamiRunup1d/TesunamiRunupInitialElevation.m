function TesunamiRunupInitialElevation
% load data
data = load('TsunamiRunupInitialCondition.mat');

% set mesh
np = 1000;
x = linspace(-500, 50000, np);
Interp = griddedInterpolant(data.x, data.eta, 'nearest');
z = Interp(x);

figure('Position',[731   352   538   295])
plot(x, z, 'r');
xlabel('x (m)', 'Interpreter', 'Latex');
ylabel('$\eta$ (m)', 'Interpreter', 'Latex');
xlim([-500, 5e4])
set(gca, 'XTick', [0:1e4:5e4], 'XTickLabel', {'0','10000', '20000', '30000', '40000', '50000'})
end% func