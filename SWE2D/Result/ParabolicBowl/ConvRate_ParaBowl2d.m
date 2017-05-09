function ConvRate_ParaBowl2d
elenum   = [80, 100, 120, 140, 160, 180];
filename = cell(numel(elenum), 1);
% triangles
order = 1;
meshtype = 'quad';
for i =1:numel(elenum)
    filename{i} = ['SWE2D_ParabolicBowl_', meshtype, '_', ...
        num2str(elenum(i)), '.nc'];
end
PostproTri  = Utilities.PostProcess.Postprocess(filename, meshtype, order);
% quadrilaterals
order = 1;
meshtype = 'quad';
for i =1:numel(elenum)
    filename{i} = ['SWE2D_ParabolicBowl_', meshtype, '_', ...
        num2str(elenum(i)), '.nc'];
end
PostproQuad = Utilities.PostProcess.Postprocess(filename, meshtype, order);

% time period
g     = 9.81;
alpha = 1.6*1e-7;
w     = sqrt(8*g*alpha);
T     = 2*pi/w;

fh = figure('color', 'w'); hold on;
markersize = 8;
linewidth = 1.5;
dofs = 4000*2./elenum';

tableTri = PostproTri.GetConvTable('h', T, @ParabolicBowl2d_ExtDepth,...
    dofs, 'M', elenum');
tableQuad = PostproQuad.GetConvTable('h', T, @ParabolicBowl2d_ExtDepth,...
    dofs, 'M', elenum');
tableTri.dofs = PostproTri.GetDofs.^2;
tableQuad.dofs = PostproQuad.GetDofs.^2;
ptri  = polyfit(log2(PostproTri.GetDofs), log2(tableTri.L2), 1)
pquad = polyfit(log2(PostproQuad.GetDofs), log2(tableQuad.L2), 1)
tableTri
tableQuad
ratioQuad = polyval(ptri, PostproQuad.GetDofs.^2)./polyval(pquad, PostproQuad.GetDofs.^2);
ratioTri = polyval(ptri, PostproTri.GetDofs.^2)./polyval(pquad, PostproTri.GetDofs.^2);
[ratioTri, ratioQuad]

plot(PostproTri.GetDofs.^2, tableTri.L2, 'bo-', 'MarkerFaceColor', 'b',...
    'MarkerSize', markersize, 'LineWidth', linewidth); hold on;
plot(PostproQuad.GetDofs.^2, tableQuad.L2, 'rs-', 'MarkerFaceColor', 'r',...
    'MarkerSize', markersize, 'LineWidth', linewidth)

tableTri = PostproTri.GetConvTable('qx', T, @ParabolicBowl2d_ExtQx,...
    dofs, 'M', elenum');
tableQuad = PostproQuad.GetConvTable('qx', T, @ParabolicBowl2d_ExtQx,...
    dofs, 'M', elenum');
tableTri
tableQuad
ptri  = polyfit(log2(PostproTri.GetDofs), log2(tableTri.L2), 1)
pquad = polyfit(log2(PostproQuad.GetDofs), log2(tableQuad.L2), 1)
% plot(PostproTri.GetDofs.^2, tableTri.L2, 'bs-', 'MarkerFaceColor', 'b',...
%     'MarkerSize', markersize, 'LineWidth', linewidth)
% plot(PostproQuad.GetDofs.^2, tableQuad.L2, 'rs-', 'MarkerFaceColor', 'r',...
%     'MarkerSize', markersize, 'LineWidth', linewidth)

tableTri = PostproTri.GetConvTable('qy', T, @ParabolicBowl2d_ExtQy,...
    dofs, 'M', elenum');
tableQuad = PostproQuad.GetConvTable('qy', T, @ParabolicBowl2d_ExtQy,...
    dofs, 'M', elenum');
tableTri
tableQuad
% plot(PostproTri.GetDofs.^2, tableTri.L2, 'b^-', 'MarkerFaceColor', 'b',...
%     'MarkerSize', markersize, 'LineWidth', linewidth)
% plot(PostproQuad.GetDofs.^2, tableQuad.L2, 'r^-', 'MarkerFaceColor', 'r',...
%     'MarkerSize', markersize, 'LineWidth', linewidth)
set(gca, 'XScale', 'log', 'YScale', 'log')
box on;
grid on;
legend({'$\rm{Tri}$', '$\rm{Quad}$'},... '$q_x$', '$q_x$', '$q_y$', '$q_y$'},...
    'Location', 'NorthEast', 'box', 'off', 'FontSize', 14,...
    'Interpreter', 'Latex')
xlabel('${DOFs}$', 'FontSize', 16, 'Interpreter', 'Latex')
ylabel('$L_2(h)$', 'FontSize', 16, 'Interpreter', 'Latex')
xlim([1.75e4, 2.e5]);
ylim([1.8e-3, 5e-3])

print -dtiff L2_DOFs.tiff -r300
end% func