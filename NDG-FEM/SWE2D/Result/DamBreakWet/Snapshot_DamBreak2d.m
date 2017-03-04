function Snapshot_DamBreak2d
%% Construct postprocess class
meshtype = 'tri';
filename = {'SWE2D_DamBreakWet_VA_tri_300.nc'};
fileID   = 1;

% create post process class for quad
Postpro = Utilities.PostProcess.Postprocess(filename, meshtype, 1);
T = 20;
time    = (0:0.005:1)*T;


figure('color', 'w');
cmap = colormap('winter');
p_h = Postpro.Snapshot2D(...
    'h',time(1)+eps, fileID,'value', ...
    cmap(end,:), cmap(1, :), [0, 10],...
    'EdgeColor', 'k',... % ������������������ɫ��͸���ȵ�
    'FaceAlpha', 0.8);
box on;
grid on;
zlim([0, 12]);
zlabel('\eta (m)','FontSize', 10);
xlabel('x (m)','FontSize', 10);
ylabel('y (m)','FontSize', 10);
view([24, 36]);
colorbar('southoutside');

x = Postpro.NcFile(fileID).GetVarData('x');
y = Postpro.NcFile(fileID).GetVarData('y');

is_Camera_on = 1; % �趨�Ƿ����ɶ���
if is_Camera_on;
    writerObj = VideoWriter([pwd,'/DamBreakWet_quad_300.avi']);
    writerObj.FrameRate=15; % �趨����֡��
    open(writerObj);	
end

for t = 1:numel(time)
    val = Postpro.GetVarData('h', time(t), fileID);
    % ���� patch ����� vertex ���ԣ��������νڵ�����˳���������
    vertex = [x(:), y(:), val(:)];
    % ͬʱ���½ڵ���ɫ
    set(p_h, 'Vertices', vertex,...
        'FaceVertexCData', val(:));
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