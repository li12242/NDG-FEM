function ComEfficient_ParaBowl2d
%% Parameter
elenum   = [20, 40, 80, 160];
fileTri  = cell(numel(elenum), 1);
fileQuad = cell(numel(elenum), 1);
degree   = 1;
for i =1:numel(elenum)
    fileTri{i}  = ['SWE2D_ParabolicBowl_', 'tri', '_', num2str(elenum(i)), '.nc'];
    fileQuad{i} = ['SWE2D_ParabolicBowl_', 'quad','_', num2str(elenum(i)), '.nc'];
end

PostproTri  = Utilities.PostProcess.Postprocess(fileTri, 'tri', degree);
PostproQuad = Utilities.PostProcess.Postprocess(fileQuad, 'quad', degree);

time  = zeros(PostproTri.nfiles, 1);
markerSize = 8;
linewidth = 1.5;

% triangle
PrintTable.dofs = PostproTri.GetDofs;
for i =1:PostproTri.nfiles
    time(i)   = PostproTri.NcFile(i).GetVarData('elapsedTime');
end% for
PrintTable.time = time*1.1;
PrintTable

figure('Color', 'w'); hold on;
plot(PrintTable.dofs.*PrintTable.dofs, PrintTable.time, 'bo-', ...
    'Markersize', markerSize, 'LineWidth', linewidth, ...
    'MarkerFaceColor', 'b')

% quadrilateral
PrintTable   = table;
PrintTable.M = elenum';
PrintTable.dofs = PostproQuad.GetDofs;
for i =1:PostproQuad.nfiles
    time(i)   = PostproQuad.NcFile(i).GetVarData('elapsedTime');
end% for
PrintTable.time = time;
PrintTable

plot(PrintTable.dofs.*PrintTable.dofs, PrintTable.time, 'rs-', ...
    'Markersize', markerSize, 'LineWidth', linewidth,...
    'MarkerFaceColor', 'r')

xlabel('$DOFs$', 'FontSize', 16, 'Interpreter', 'Latex')
ylabel('$\rm{Elapsed \: Time (s)}$', 'FontSize', 16, 'Interpreter', 'Latex')

box on; grid on;
legend({'$\rm{Tri}$', '$\rm{Quad}$'}, 'box', 'off', 'FontSize', 16,...
    'Location', 'NorthWest', 'Interpreter', 'Latex');
end% func