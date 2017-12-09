function Animat_TsunamiRunup

%% Parameter
T        = 22;
meshtype = 'quad';
filepath = pwd;
filename = {[filepath,'/SWE2D_TsuamiRunup_quad_200_5e4.nc']};
order    = 1;
Postpro  = Utilities.PostProcess.Postprocess(filename, meshtype, order);
fileID   = 1;
time     = linspace(1e-12, T, 100);

x = Postpro.NcFile(fileID).GetVarData('x');
y = Postpro.NcFile(fileID).GetVarData('y');

figure('color', 'w'); hold on;
bh  = Postpro.SnapshotConst2D('bot',fileID);
set(bh, 'FaceColor', [0.8, 0.8, .8]);
bot = Postpro.NcFile(1).GetVarData('bot');

cmap = colormap('winter');
ph = Postpro.Snapshot2D(...
    'h',time(1), 1,'value', ...
    cmap(end,:), cmap(1, :), [-0., 0.1],...
    'EdgeColor', 'none',... % 其他参数，包括边颜色，透明度等
    'FaceAlpha', 0.8);
zlim([-0.15, .15]);
view(-140, 50);
zlabel('elevation (m)','FontSize', 14);
xlabel('x (m)','FontSize', 16);
ylabel('y (m)','FontSize', 16);
box on;
grid on;
colorbar('southoutside');

is_Camera_on = 1; % 设定是否生成动画
if is_Camera_on;
    writerObj = VideoWriter([pwd,'/TsunamiRunup_quad_200.avi']);
    writerObj.FrameRate=15; % 设定动画帧率
    open(writerObj);	
end

for i = 1:numel(time)
    
    val = Postpro.GetVarData('h', time(i), fileID);
    val(val<1e-3) = nan;
    % 更新 patch 对象的 vertex 属性，否则多边形节点连接顺序产生错误
    vertex = [x(:), y(:), val(:)+bot(:)];
    % 同时更新节点颜色
    set(ph, 'Vertices', vertex,...
        'FaceVertexCData', val(:)+bot(:));
    
    if is_Camera_on
        frame = getframe(gcf);
        writeVideo(writerObj,frame); 
    end% if
end% for

if is_Camera_on
    close(writerObj);
end
end% func