function Convection2d_Snapshot
%CONVECTION2D_SNAPSHOT Snapshot of 2-d convection problem

%% Parameter
T        = 2.4;
meshtype = 'quad';
order    = 1;
filename = {'Convection2D_quad_1_100.nc'};
PostproTri  = Utilities.PostProcess.Postprocess(filename, meshtype, order);
fileID   = 1;
time     = [(1e-12:0.1:1)*T, T];

%% draw pic
for i = 1:numel(time)
    figure
    PostproTri.Snapshot2D('var', time(i), fileID);
    zlim([-.5, 1.2]);
    view(30, 32);
    zlabel('var','FontSize', 14);
    xlabel('x'  ,'FontSize', 14);
    ylabel('y'  ,'FontSize', 14);
    box on
end% for

% time = PostproTri.NcFile(1).GetVarData('time');
% for i = 1:5:numel(time)
%     var = PostproTri.GetVarData('var', time(i), fileID);
%     fprintf('time=%f, max var=%f\n',time(i), max(max(var)));
% end
end

