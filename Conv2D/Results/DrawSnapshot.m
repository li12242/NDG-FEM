function DrawSnapshot
%% parameters
T        = 0.5;
nt = 100;
time     = linspace(eps, T, nt);
fileID   = 1; % index of result file to draw
ele      = 50;
deg      = 2;
meshtype = 'quad';
filename = cell(numel(ele), 1);
for i = 1:numel(ele)
    filename{i} = ['Convection2D_', meshtype, '_', ...
        num2str(deg),'_',num2str(ele(i)), '.nc'];
end% for

%% construct postprocess class
Postpro = Utilities.PostProcess.Postprocess(filename, meshtype, deg);

%% draw snapshots
x = Postpro.NcFile(fileID).GetVarData('x');
y = Postpro.NcFile(fileID).GetVarData('y');

p_h = Postpro.Snapshot2D('p', time(1), fileID, 'default');
view(-29, 20);
xlim([-1, 1]); ylim([-1, 1]); 
%zlim([-0.1, 1.5]);
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
zlabel('$C$', 'Interpreter', 'latex')

for i = 1:numel(time)
    val = Postpro.GetVarData('p', time(i), fileID);
    
    vertex = [x(:), y(:), val(:)];
    % 同时更新节点颜色
    set(p_h, 'Vertices', vertex,...
        'FaceVertexCData', val(:));
    drawnow;
    
end
end% func