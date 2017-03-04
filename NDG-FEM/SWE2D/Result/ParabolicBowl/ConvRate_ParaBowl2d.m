function ConvRate_ParaBowl2d
elenum   = [81, 100, 120, 140, 160, 180];
filename = cell(numel(elenum), 1);
% triangles
order = 1;
meshtype = 'quad';
for i =1:numel(elenum)
    filename{i} = ['swe2d_', meshtype, '_', num2str(order) ,'_', ...
        num2str(elenum(i)), '.0-1.nc'];
end
PostproTri  = Utilities.PostProcess.Postprocess(filename, meshtype, order);
% quadrilaterals
order = 1;
meshtype = 'quad';
for i =1:numel(elenum)
    filename{i} = ['swe2d_', meshtype, '_', num2str(order) ,'_', ...
        num2str(elenum(i)), '.0-1.nc'];
end
PostproQuad = Utilities.PostProcess.Postprocess(filename, meshtype, order);

% time period
g     = 9.81;
alpha = 1.6*1e-7;
w     = sqrt(8*g*alpha);
T     = 2*pi/w;

fh = figure('color', 'w'); hold on;
box on;
grid on;
markersize = 6;
linewidth = 1.5;
dofs = 4000*2./elenum';

tableTri = PostproTri.GetConvTable('h', T, @ParabolicBowl2d_ExtDepth,...
    dofs, 'M', elenum');
tableQuad = PostproQuad.GetConvTable('h', T, @ParabolicBowl2d_ExtDepth,...
    dofs, 'M', elenum');
tableTri
tableQuad
plot(PostproTri.GetDofs.^2, tableTri.L2, 'bo-', 'MarkerFaceColor', 'b',...
    'MarkerSize', markersize, 'LineWidth', linewidth); hold on;
plot(PostproQuad.GetDofs.^2, tableQuad.L2, 'ro-', 'MarkerFaceColor', 'r',...
    'MarkerSize', markersize, 'LineWidth', linewidth)

tableTri = PostproTri.GetConvTable('qx', T, @ParabolicBowl2d_ExtQx,...
    dofs, 'M', elenum');
tableQuad = PostproQuad.GetConvTable('qx', T, @ParabolicBowl2d_ExtQx,...
    dofs, 'M', elenum');
tableTri
tableQuad
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
ah = fh.Children;
set(ah, 'XScale', 'log', 'YScale', 'log')
legend({'$\rm{Tri}$', '$\rm{Quad}$'},... '$q_x$', '$q_x$', '$q_y$', '$q_y$'},...
    'Location', 'NorthEast', 'box', 'off', 'FontSize', 14,...
    'Interpreter', 'Latex')
xlabel('${DOFs}$', 'FontSize', 16, 'Interpreter', 'Latex')
ylabel('$L_2(h)$', 'FontSize', 16, 'Interpreter', 'Latex')
xlim([1e4,1e6]);
ylim([0.015, 0.045])

print -dtiff L2_DOFs.tiff -r300
end% func