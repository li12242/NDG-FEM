function Snapshot_ParaBowl2d
%SNAPSHOT_PARABOWL2D Summary of this function goes here
%   Detailed explanation goes here
% Parameter
T        = 20;
% casename = 'SWE2D_FlowOver3BumpsUniform_';
meshtype = 'quad';
% filename = {[casename, meshtype, '_100.nc']};
filename = {'SWE2D_ParabolicBowl_quad_80.nc'};
order    = 1;
Postpro  = Utilities.PostProcess.Postprocess(filename, meshtype, order);
fileID   = 1;
time     = (0:0.005:1)*T;

figure('color', 'w')
% plot bottom topography
bot = Postpro.NcFile(1).GetVarData('bot');
bh = Postpro.SnapshotConst2D('bot',fileID);
set(bh, 'FaceColor', [.8, .8, .8]);

hold on;
box on;
grid on;

% plot water level
cmap = colormap('winter');
p_h = Postpro.Snapshot2D(...
    'h',time(1)+eps, 1,'value', ...
    cmap(end,:), cmap(1, :), [0, 1.875],...
    'EdgeColor', 'none',... % 其他参数，包括边颜色，透明度等
    'FaceAlpha', 0.8);

% set axes properties
set(gca, 'DataAspectRatio', [1000,1000,1])
xlim([0, 4000])
zlim([0, 5]);
zlabel('\eta (m)','FontSize', 10);
xlabel('x (m)','FontSize', 10);
ylabel('y (m)','FontSize', 10);
view(-40, 24);

is_Camera_on = 1; % 设定是否生成动画
if is_Camera_on;
    writerObj = VideoWriter([pwd,'/ParaBowl2d_quad_80.avi']);
    writerObj.FrameRate=15; % 设定动画帧率
    open(writerObj);	
end

% Update water level
x = Postpro.NcFile(fileID).GetVarData('x');
y = Postpro.NcFile(fileID).GetVarData('y');
for i = 1:numel(time)
    val = Postpro.GetVarData('h', time(i), fileID);
    val(val<1e-2) = nan;
    % 更新 patch 对象的 vertex 属性，否则多边形节点连接顺序产生错误
    vertex = [x(:), y(:), val(:)+bot(:)];
    % 同时更新节点颜色
    set(p_h, 'Vertices', vertex,...
        'FaceVertexCData', val(:)+bot(:));
    drawnow;
    
    if is_Camera_on
        frame = getframe(gcf);
        writeVideo(writerObj,frame); 
    end% if
end% for

if is_Camera_on
    close(writerObj);
end

end

