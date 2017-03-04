function SnapShot_TsunamiRunup

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
    'EdgeColor', 'none',... % ������������������ɫ��͸���ȵ�
    'FaceAlpha', 0.8);
zlim([-0.15, .15]);
view(-140, 50);
zlabel('elevation (m)','FontSize', 14);
xlabel('x (m)','FontSize', 16);
ylabel('y (m)','FontSize', 16);
box on;
grid on;
colorbar('southoutside');

is_Camera_on = 1; % �趨�Ƿ����ɶ���
if is_Camera_on;
    writerObj = VideoWriter([pwd,'/TsunamiRunup_quad_200.avi']);
    writerObj.FrameRate=15; % �趨����֡��
    open(writerObj);	
end

for i = 1:numel(time)
    
    val = Postpro.GetVarData('h', time(i), fileID);
    val(val<1e-3) = nan;
    % ���� patch ����� vertex ���ԣ��������νڵ�����˳���������
    vertex = [x(:), y(:), val(:)+bot(:)];
    % ͬʱ���½ڵ���ɫ
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