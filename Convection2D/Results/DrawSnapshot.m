function DrawSnapshot
%% parameters
T        = 2.4;
time     = [(1e-12:1/4:1)*T, T];
fileID   = 1; % index of result file to draw
ele      = [100];
deg      = 1;
meshtype = 'quad';
filename = cell(numel(ele), 1);
for i = 1:numel(ele)
    filename{i} = ['Convection2D_', meshtype, '_', ...
        num2str(deg),'_',num2str(ele(i)), '.nc'];
end% for

%% construct postprocess class
PostproConv2d = Utilities.PostProcess.Postprocess(filename, meshtype, deg);

%% draw snapshots
for i = 1:numel(time)
    figure;
    PostproConv2d.Snapshot2D('var', time(i), fileID);
    view(-29, 20)
    title(['$t=',num2str(time(i)),'$'], 'Interpreter', 'latex',...
        'FontSize', 18)
    
    xlim([0, 1]); ylim([0, 1]); zlim([-0.1, 1.5]);
    xlabel('$x$', 'Interpreter', 'latex')
    ylabel('$y$', 'Interpreter', 'latex')
    zlabel('$C$', 'Interpreter', 'latex')
end
end% func