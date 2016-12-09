function Snapshot_DamBreak2d
%% Construct postprocess class
meshtype = 'tri';
filename = {'SWE2D_DamBreakWet_VA_tri_300.nc'};
fileID   = 1;

% create post process class for quad
Postpro = Utilities.PostProcess.Postprocess(filename, meshtype, 1);

time = Postpro.NcFile(1).GetVarData('time');

cmap = colormap('winter');
p_h = Postpro.Snapshot2D(...
    'h',time(1)+eps, 1,'value', ...
    cmap(end,:), cmap(1, :), [0, 10],...
    'EdgeColor', 'k',... % 其他参数，包括边颜色，透明度等
    'FaceAlpha', 0.8);
box on;
grid on;
view([24, 36]);

x = Postpro.NcFile(fileID).GetVarData('x');
y = Postpro.NcFile(fileID).GetVarData('y');
for t = 1:numel(time)
    val = Postpro.NcFile(1).GetTimeVarData('h', t);
    % 更新 patch 对象的 vertex 属性，否则多边形节点连接顺序产生错误
    vertex = [x(:), y(:), val(:)];
    % 同时更新节点颜色
    set(p_h, 'Vertices', vertex,...
        'FaceVertexCData', val(:));
    drawnow;
end% for
end