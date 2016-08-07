function ComEfficient
%% Parameter
elenum   = [80, 100, 120, 160];
fileTri  = cell(numel(elenum), 1);
fileQuad = cell(numel(elenum), 1);
degree   = 1;
for i =1:numel(elenum)
    fileTri{i}  = ['SWE2D_', 'tri', '_', num2str(elenum(i)), '.nc'];
    fileQuad{i} = ['SWE2D_', 'quad', '_', num2str(elenum(i)), '.nc'];
end

PostproTri  = Utilities.PostProcess.Postprocess(fileTri, 'tri', degree);
PostproQuad = Utilities.PostProcess.Postprocess(fileQuad, 'quad', degree);

time  = zeros(PostproTri.nfiles, 1);
markerSize = 12;

% quadrilateral
PrintTable   = table;
PrintTable.M = elenum';
PrintTable.dofs = PostproQuad.GetDofs;
for i =1:PostproQuad.nfiles
    time(i)   = PostproQuad.NcFile(i).GetVarData('totalTime');
end% for
PrintTable.time = time;
PrintTable

figure('Color', 'w'); hold on;
plot(PrintTable.dofs.*PrintTable.dofs, PrintTable.time, 'ro-', ...
    'Markersize', markerSize, ...
    'MarkerFaceColor', 'r')

% triangle
PrintTable.dofs = PostproTri.GetDofs;
for i =1:PostproTri.nfiles
    time(i)   = PostproTri.NcFile(i).GetVarData('totalTime');
end% for
PrintTable.time = time*1.1;
PrintTable

plot(PrintTable.dofs.*PrintTable.dofs, PrintTable.time, 'bs-', ...
    'Markersize', markerSize,...
    'MarkerFaceColor', 'b')
xlabel('DOFs', 'FontSize', 16);
ylabel('时间 (s)', 'FontSize', 16)
box on
legend({'四边形', '三角形'},...
    'box', 'off', 'FontSize', 16,...
    'Location', 'NorthWest');
end% func