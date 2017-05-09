function ComEfficient_ParaBowl2d
%% Parameter
elenum   = [80, 100, 120, 140];
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
PrintTableTri   = table;
PrintTableTri.M = elenum';
PrintTableTri.dofs = PostproTri.GetDofs;
for i =1:PostproTri.nfiles
    time(i) = PostproTri.NcFile(i).GetVarData('elapse');
end% for
PrintTableTri.time = time;
PrintTableTri

figure('Color', 'w'); hold on;
plot(PrintTableTri.dofs.*PrintTableTri.dofs, PrintTableTri.time, 'bo-', ...
    'Markersize', markerSize, 'LineWidth', linewidth, ...
    'MarkerFaceColor', 'b')

% quadrilateral
PrintTableQuad   = table;
PrintTableQuad.M = elenum';
PrintTableQuad.dofs = PostproQuad.GetDofs;
for i =1:PostproQuad.nfiles
    time(i) = PostproQuad.NcFile(i).GetVarData('elapse');
end% for
PrintTableQuad.time = time;
PrintTableQuad

ptri  = polyfit(PrintTableTri.dofs.^2, PrintTableTri.time, 1);
pquad = polyfit(PrintTableQuad.dofs.^2, PrintTableQuad.time, 1);
ratioQuad = polyval(ptri, PrintTableQuad.dofs.^2)./polyval(pquad, PrintTableQuad.dofs.^2);
ratioTri = polyval(ptri, PrintTableTri.dofs.^2)./polyval(pquad, PrintTableTri.dofs.^2);
[ratioTri, ratioQuad]

plot(PrintTableQuad.dofs.*PrintTableQuad.dofs, PrintTableQuad.time, 'rs-', ...
    'Markersize', markerSize, 'LineWidth', linewidth,...
    'MarkerFaceColor', 'r')
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('$DOFs$', 'FontSize', 16, 'Interpreter', 'Latex')
ylabel('$\rm{Elapsed \,\, Time (s)}$', 'FontSize', 16, 'Interpreter', 'Latex')
xlim([1.75e4, 2.e5]);
% ylim([2e1, 6e2]);

box on; grid on;
legend({'$\rm{Tri}$', '$\rm{Quad}$'}, 'box', 'off', 'FontSize', 16,...
    'Location', 'NorthWest', 'Interpreter', 'Latex');
print -dtiff L2_DOFs.tiff -r300
end% func