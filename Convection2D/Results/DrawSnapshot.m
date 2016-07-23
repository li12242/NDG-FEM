function DrawSnapshot
%% parameters
T        = 2.4;
time     = (0:1/4:1)*T;
fileID   = 2; % index of result file to draw
ele      = [20, 40, 60, 80];
deg      = 1;
meshtype = 'quad';
filename = cell(numel(ele), 1);
for i = 1:numel(ele)
    filename{i} = ['Convection2D_', meshtype, '_', num2str(deg),'_',num2str(ele(i)), '.nc'];
end% for

%% construct postprocess class
PostproConv2d = Utilities.PostProcess.Postprocess(filename, meshtype, deg);

%% draw snapshots
for i = 1:numel(time)
    figure;
    PostproConv2d.Snapshot2D('var', time(i), fileID)
    view(-29, 20)
    title(['$t=',num2str(time(i)),'$'], 'Interpreter', 'latex',...
        'FontSize', 18)
    
    xlim([-1, 1]); ylim([-1, 1]); zlim([-0.2, 1.]);
    xlabel('$x$', 'Interpreter', 'latex')
    ylabel('$y$', 'Interpreter', 'latex')
    zlabel('$C$', 'Interpreter', 'latex')
end
end% func