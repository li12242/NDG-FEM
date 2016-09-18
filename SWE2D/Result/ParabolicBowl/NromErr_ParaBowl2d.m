function ParaBowl2d_NromErr
%% parameters
meshtype = 'quad';
elenum   = [80, 100, 120, 160, 200];
filename = cell(numel(elenum), 1);
for i =1:numel(elenum)
    filename{i} = ['SWE2D_', meshtype, '_', num2str(elenum(i)), '.nc'];
end

% create post process class for quad
PostproQuad = Utilities.PostProcess.Postprocess(filename, meshtype, 1);

% parameters
meshtype = 'tri';
for i =1:numel(elenum)
    filename{i} = ['SWE2D_', meshtype, '_', num2str(elenum(i)), '.nc'];
end

% create post process class for quad
PostproTri = Utilities.PostProcess.Postprocess(filename, meshtype, 1);

% time
g     = 9.81;
alpha = 1.6*1e-7;
w     = sqrt(8*g*alpha);
T     = 2*pi/w;

%% plot
figure; ah1 = axes; hold on
xlabel('$\sqrt{DOFs}$', 'Interpreter', 'Latex');
ylabel('$L_2$', 'Interpreter', 'Latex')
figure; ah2 = axes; hold on
xlabel('$\sqrt{DOFs}$', 'Interpreter', 'Latex');
ylabel('$L_{\infty}$', 'Interpreter', 'Latex');
%% calculate error for quad
% error for H
PrintTable   = table;
PrintTable.M = elenum';
PrintTable.dofs = PostproQuad.GetDofs;

[errL2, errLinf]= GetErr(PostproQuad, 'h', T, @ParabolicBowl2d_ExtDepth);
PrintTable.L2   = errL2;
PrintTable.Linf = errLinf;
fprintf('Quad - water depth');
PrintTable
[p1, p2] = DrawErr(ah1, ah2, PrintTable, 'b','o');

% error for qx
[errL2, errLinf]= GetErr(PostproQuad, 'qx', T, @ParabolicBowl2d_ExtQx);
PrintTable.L2   = errL2;
PrintTable.Linf = errLinf;
fprintf('Quad - flow flux x');
PrintTable
[p3, p4] = DrawErr(ah1, ah2, PrintTable, 'b','+');

%% get error for triangle
PrintTable.dofs = PostproTri.GetDofs;
[errL2, errLinf]= GetErr(PostproTri, 'h', T, @ParabolicBowl2d_ExtDepth);
PrintTable.L2   = errL2;
PrintTable.Linf = errLinf;
fprintf('Quad - water depth');
PrintTable
[p5, p6] = DrawErr(ah1, ah2, PrintTable, 'r','s');

% error for qx
[errL2, errLinf]= GetErr(PostproTri, 'qx', T, @ParabolicBowl2d_ExtQx);
PrintTable.L2   = errL2;
PrintTable.Linf = errLinf;
fprintf('Quad - flow flux x');
PrintTable
[p7, p8] = DrawErr(ah1, ah2, PrintTable, 'r','^');


%% set properties
set(ah1,'XScale','log','YScale','log',...
        'XGrid' ,'on', 'YGrid', 'on',...
        'Box','on',...
        'XLim', [150, 500], 'YLim', [1.5e-3, 2e-2]);
set(ah2,'XScale','log','YScale','log',...
        'XGrid' ,'on', 'YGrid', 'on',...
        'Box','on',...
        'XLim', [150, 500], 'YLim', [2e-2, 0.1]);
% legend
legend(ah1, [p1, p3, p5, p7],....
    {'Quad - H', 'Quad - Qx', 'Tri - H', 'Tri - Qx'}, ...
    'box', 'off', 'Interpreter', 'Latex');
legend(ah2, [p2, p4, p6, p8],...
    {'Quad - H', 'Quad - Qx', 'Tri - H', 'Tri - Qx'}, ....
    'box', 'off', 'Interpreter', 'Latex');
end% func

function [p1, p2] = DrawErr(ah1, ah2, PrintTable, Color, Symbol)
p1 = plot(ah1, PrintTable.dofs, PrintTable.L2, [Color,Symbol,'-'],...
    'Markersize', 8, 'MarkerFaceColor', Color);
p2 = plot(ah2, PrintTable.dofs, PrintTable.Linf, [Color,Symbol,'-'],...
    'Markersize', 8, 'MarkerFaceColor', Color);
end% func

function [errL2, errLinf] = GetErr(PostproSWE2D, varname, T, exFunH)
errL2        = zeros(PostproSWE2D.nfiles, 1);
errLinf      = zeros(PostproSWE2D.nfiles, 1);

for i =1:PostproSWE2D.nfiles
    errL2(i)   = PostproSWE2D.NormErr(varname, T, exFunH, 'L2', i);
    errLinf(i) = PostproSWE2D.NormErr(varname, T, exFunH, 'Linf', i);
end% for
end% func
