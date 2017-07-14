function LakeAtRestVertify
filename = 'SWE1D_LakeAtRest.nc';
filename1 = 'SWE1D_LakeAtRest1.nc'; % h-type refinement
time = ncread(filename, 'time');
time1 = ncread(filename1, 'time');
% timestep = numel(time)
x = ncread(filename, 'x');

% bottom topography
bedElevation = zeros(size(x));
a = 1.2; rm = 0.4; r = abs(x - 0.5);
index = (r < rm);
bedElevation(index) = a*exp(-0.5./(rm.^2 - r(index).^2))./exp(-0.5./rm^2);

FinalTime = [0.000, 0.2];

for itime = 1:numel(FinalTime)
    [~, index] = min(abs(time - FinalTime(itime)));
    
    h = ncread(filename, 'h', [1, index],[inf, 1]);
    q = ncread(filename, 'q', [1, index],[inf, 1]);
    
%     [~, index1] = min(abs(time1 - FinalTime(itime)));
    index1 = 1;
    
    h1 = ncread(filename1, 'h', [1, index1],[inf, 1]);
    q1 = ncread(filename1, 'q', [1, index1],[inf, 1]);
    
    % draw picture
    figure('Position', [1022 637 524 193]);
    plot(x, h+bedElevation, 'b-.', 'LineWidth', 2); hold on;
    plot(x, h1+bedElevation, 'r--', 'LineWidth', 1.5)
    plot(x, bedElevation, 'k')
    xlim([0.384, 0.426]); ylim([0.93,1.1])
    xlabel('x', 'Interpreter', 'Latex');
    ylabel('$\eta$', 'Interpreter', 'Latex');

    figure('Position', [1022 637 524 193]);
    plot(x, q, 'b'); hold on;
    plot(x, q1, 'r');
    xlabel('x', 'Interpreter', 'Latex');
    ylabel('q', 'Interpreter', 'Latex');
end% for
figure(1);
t = legend('CONV', 'HREF');
set(t, 'box', 'off');
end% func