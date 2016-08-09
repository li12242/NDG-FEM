function DamBreakDry1d_Snapshot
%DAMBREAK1D_SNAPSHOT Snapshots of one dimensional DamBreak problem
%   The numerical solution is compared with the exact solution at various
%   time.

% parameter
T = 20;
meshType = 'line';
casename = 'DamBreakDry';
nele     = [100, 200, 300, 400];
filename = cell(numel(nele), 1);
for i = 1:numel(nele)
    filename{i} = ['SWE1D_', casename, '_', num2str(nele(i)),'.nc'];
end
Postpro  = Utilities.PostProcess.Postprocess(filename, meshType, 1);
fileID   = 4;
time     = [4,8,12];

ne       = 400; % number of exact solution points
xe       = linspace(0, 1e3, ne);

% draw results
figure('Position', [519   578   534   189]); hold on;
colorRBG  = {[0, 0.4000, 1.0000], [1, 0, 0], [0., 0.7980, 0]};
labelSize = 15;
legendSize = 15;
for i = 1:numel(time)
    p{i} = Postpro.Snapshot1D('h', time(i), fileID,...
        'Marker','o', 'Markersize', 5, ...
        'Color', colorRBG{i});
    he   = DamBreakDry1d_ExactH(xe, time(i));
    pe   = plot(xe, he, 'k-');
    ylim([0, 10]);
    xlabel('$\mathrm{x} \,(m)$',...
        'FontSize', labelSize,...
        'Interpreter', 'latex');
    ylabel('$\mathrm{Depth} \,(m)$',...
        'FontSize', labelSize, ...
        'Interpreter', 'latex');
    box on;
%     grid on;
end% for
t = legend([p{1}(1), p{2}(1), p{3}(1), pe], ...
    't=4s', 't=8s', 't=12s', 'Exact');
t.Box = 'off';
t.FontSize = legendSize;


figure('Position', [519   578   534   189]); hold on;
for i = 1:numel(time)
    p{i} = Postpro.Snapshot1D('q', time(i), fileID,...
        'Marker','o', 'Markersize', 5, ...
        'Color', colorRBG{i});
    he   = DamBreakDry1d_ExactQ(xe, time(i));
    pe   = plot(xe, he, 'k-');
    ylim([0, 30]);
    xlabel('$\mathrm{x} \,(m)$',...
        'FontSize', labelSize, ...
        'Interpreter', 'latex');
    ylabel('$\mathrm{Flux} \,(m^2/s)$',...
        'FontSize', labelSize, ...
        'Interpreter', 'latex');
    box on;
%     grid on;
end% for
t = legend([p{1}(1), p{2}(1), p{3}(1), pe], ...
    {'t=4s', 't=8s', 't=12s', 'Exact'});
t.Box = 'off';
t.FontSize = legendSize;
end% func