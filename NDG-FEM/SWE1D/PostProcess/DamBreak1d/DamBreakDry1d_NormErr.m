function DamBreakDry1d_NormErr
%DAMBREAK1D_NORMERR Calculation of the L2 and Linf error.
%   Compare the numerical solution with the exact solution in L2 and Linf
%   norm error, and plot the error over the number of DOFs.

casename = 'DamBreakDry';
nele     = [100, 200, 300, 400];
filename = cell(numel(nele), 1);
for i = 1:numel(nele)
    filename{i} = ['SWE1D_', casename, '_', num2str(nele(i)),'.nc'];
end
T        = 20;
order    = 1;
PostPro  = Utilities.PostProcess.Postprocess(filename, 'line', order);

%% plot
figure; ah1 = axes; hold on
xlabel('$\sqrt{DOFs}$', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('$L_2$',         'Interpreter', 'Latex', 'FontSize', 15)
figure; ah2 = axes; hold on
xlabel('$\sqrt{DOFs}$', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('$L_{\infty}$',  'Interpreter', 'Latex', 'FontSize', 15);


PrintTable   = table;
PrintTable.M = nele';
PrintTable.dofs = PostPro.GetDofs;

[errL2, errLinf] = GetErr(PostPro, 'h', T, @DamBreakDry1d_ExactH);
PrintTable.L2   = errL2;
PrintTable.Linf = errLinf;
fprintf('water depth');
PrintTable
[p1, p2] = DrawErr(ah1, ah2, PrintTable, 'b', 'o');

[errL2, errLinf] = GetErr(PostPro, 'q', T, @DamBreakDry1d_ExactQ);
PrintTable.L2   = errL2;
PrintTable.Linf = errLinf;
fprintf('flux');
PrintTable
[p3, p4] = DrawErr(ah1, ah2, PrintTable, 'r', 'o');

set(ah1,'XScale','log','YScale','log','XGrid' ,'on',...
    'YGrid', 'on','Box','on','XLim', [150, 1e3])
set(ah2,'XScale','log','YScale','log','XGrid' ,'on',...
    'YGrid', 'on','Box','on','XLim', [150, 1e3])

legend(ah1, [p1, p3], {'h', 'q'}, ...
    'Interpreter', 'Latex', 'Box', 'off', 'FontSize', 15);
legend(ah2, [p2, p4], {'h', 'q'}, ...
    'Interpreter', 'Latex', 'Box', 'off', 'FontSize', 15);

end

function [p1, p2] = DrawErr(ah1, ah2, PrintTable, Color, Symbol)
p1 = plot(ah1, PrintTable.dofs, PrintTable.L2, [Color,Symbol,'-'],...
    'Markersize', 8, 'MarkerFaceColor', Color);
p2 = plot(ah2, PrintTable.dofs, PrintTable.Linf, [Color,Symbol,'-'],...
    'Markersize', 8, 'MarkerFaceColor', Color);
end% func

function [errL2, errLinf] = GetErr(PostproSWE1D, varname, T, exFunH)
errL2        = zeros(PostproSWE1D.nfiles, 1);
errLinf      = zeros(PostproSWE1D.nfiles, 1);

for i =1:PostproSWE1D.nfiles
    errL2(i)   = PostproSWE1D.NormErr(varname, T, exFunH, 'L2', i);
    errLinf(i) = PostproSWE1D.NormErr(varname, T, exFunH, 'Linf', i);
end% for
end

