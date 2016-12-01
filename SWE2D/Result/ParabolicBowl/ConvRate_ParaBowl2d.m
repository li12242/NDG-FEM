function ConvRate_ParaBowl2d
elenum   = [20, 40, 80, 160];
filename = cell(numel(elenum), 1);
% triangles
meshtype = 'tri';
for i =1:numel(elenum)
    filename{i} = ['SWE2D_ParabolicBowl_', meshtype, '_', ...
        num2str(elenum(i)), '.nc'];
end
PostproTri  = Utilities.PostProcess.Postprocess(filename, 'tri', 1);
% quadrilaterals
meshtype = 'quad';
for i =1:numel(elenum)
    filename{i} = ['SWE2D_ParabolicBowl_', meshtype, '_', ...
        num2str(elenum(i)), '.nc'];
end
PostproQuad = Utilities.PostProcess.Postprocess(filename, 'quad', 1);

% time period
g     = 9.81;
alpha = 1.6*1e-7;
w     = sqrt(8*g*alpha);
T     = 2*pi/w;

figure('color', 'w'); hold on;
box on;
grid on;
markersize = 6;
linewidth = 1.5;

tableTri = PostproTri.GetConvTable('h', T, @ParabolicBowl2d_ExtDepth,...
    'M', elenum');
tableQuad = PostproQuad.GetConvTable('h', T, @ParabolicBowl2d_ExtDepth,...
    'M', elenum');
tableTri
plot(tableTri.dofs, tableTri.L2, 'bo-', 'MarkerFaceColor', 'b',...
    'MarkerSize', markersize, 'LineWidth', linewidth)
plot(tableQuad.dofs, tableQuad.L2, 'ro-', 'MarkerFaceColor', 'r',...
    'MarkerSize', markersize, 'LineWidth', linewidth)

tableTri = PostproTri.GetConvTable('qx', T, @ParabolicBowl2d_ExtQx,...
    'M', elenum');
tableQuad = PostproQuad.GetConvTable('qx', T, @ParabolicBowl2d_ExtQx,...
    'M', elenum');
tableTri
plot(tableTri.dofs, tableTri.L2, 'bs-', 'MarkerFaceColor', 'b',...
    'MarkerSize', markersize, 'LineWidth', linewidth)
plot(tableQuad.dofs, tableQuad.L2, 'rs-', 'MarkerFaceColor', 'r',...
    'MarkerSize', markersize, 'LineWidth', linewidth)

tableTri = PostproTri.GetConvTable('qy', T, @ParabolicBowl2d_ExtQy,...
    'M', elenum');
tableQuad = PostproQuad.GetConvTable('qy', T, @ParabolicBowl2d_ExtQy,...
    'M', elenum');
tableTri
plot(tableTri.dofs, tableTri.L2, 'b^-', 'MarkerFaceColor', 'b',...
    'MarkerSize', markersize, 'LineWidth', linewidth)
plot(tableQuad.dofs, tableQuad.L2, 'r^-', 'MarkerFaceColor', 'r',...
    'MarkerSize', markersize, 'LineWidth', linewidth)

set(gca, 'XScale', 'log', 'YScale', 'log')
legend({'$h (\rm{Quad})$', '$h (\rm{Tri})$', '$q_x$', '$q_x$', '$q_y$', '$q_y$'},...
    'Location', 'NorthWest', 'box', 'off', 'FontSize', 14,...
    'Interpreter', 'Latex')
xlabel('$\sqrt{DOFs}$', 'FontSize', 16, 'Interpreter', 'Latex')
ylabel('$L_2$', 'FontSize', 16, 'Interpreter', 'Latex')

end% func