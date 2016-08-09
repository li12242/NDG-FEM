function DamBreakDry1d_ConvRate
%DAMBREAK1D_CONVRATE Computation of the convergence rate.
%   Detailed explanation goes here

casename = 'DamBreakDry';
nele     = [100, 200, 300, 400];
filename = cell(numel(nele), 1);
for i = 1:numel(nele)
    filename{i} = ['SWE1D_', casename, '_', num2str(nele(i)),'.nc'];
end
T        = 20;
order    = 1;
PostProSWE1d  = Utilities.PostProcess.Postprocess(filename, 'line', order);

PrintTable   = table;
PrintTable.M = nele';
PrintTable.dofs = PostProSWE1d.GetDofs;

[errL2, errLinf] = GetErr(PostProSWE1d, 'h', T, @DamBreakDry1d_ExactH);
a2   = PostProSWE1d.ConvRate('h', T, @DamBreakDry1d_ExactH, 'L2');
ainf = PostProSWE1d.ConvRate('h', T, @DamBreakDry1d_ExactH, 'Linf');

PrintTable.L2   = errL2;
PrintTable.a2   = a2;
PrintTable.Linf = errLinf;
PrintTable.ainf = ainf;
fprintf('water depth');
PrintTable

[errL2, errLinf] = GetErr(PostProSWE1d, 'q', T, @DamBreakDry1d_ExactQ);
a2   = PostProSWE1d.ConvRate('q', T, @DamBreakDry1d_ExactQ, 'L2');
ainf = PostProSWE1d.ConvRate('q', T, @DamBreakDry1d_ExactQ, 'Linf');

PrintTable.L2   = errL2;
PrintTable.a2   = a2;
PrintTable.Linf = errLinf;
PrintTable.ainf = ainf;
fprintf('flux');
PrintTable

end

function [errL2, errLinf] = GetErr(PostproSWE1D, varname, T, exFunH)
errL2        = zeros(PostproSWE1D.nfiles, 1);
errLinf      = zeros(PostproSWE1D.nfiles, 1);

for i =1:PostproSWE1D.nfiles
    errL2(i)   = PostproSWE1D.NormErr(varname, T, exFunH, 'L2', i);
    errLinf(i) = PostproSWE1D.NormErr(varname, T, exFunH, 'Linf', i);
end% for
end% func
