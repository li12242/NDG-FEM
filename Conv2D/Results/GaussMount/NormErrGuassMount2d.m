function NormErrGuassMount2d
%NORMERRGUASSMOUNT2D Summary of this function goes here
%   Detailed explanation goes here

%% parameters
meshtype = 'quad';
N        = 1;
elenum   = [50];
filename = cell(numel(elenum), 1);
for i =1:numel(elenum)
    filename{i} = ['Convection2D_', ...
        meshtype, '_', num2str(N),'_',num2str(elenum(i)),'.nc'];
end

Postpro = Utilities.PostProcess.Postprocess(filename, meshtype, N);
T = 2.4;
exFunH = @ExactGaussMount2d;

%% table for norm error
errL2        = zeros(Postpro.nfiles, 1);
errLinf      = zeros(Postpro.nfiles, 1);

PrintTable   = table;
PrintTable.M = elenum';
PrintTable.dofs = Postpro.GetDofs;

for i =1:Postpro.nfiles
    errL2(i)   = Postpro.NormErr('var', T, exFunH, 'L2', i);
    errLinf(i) = Postpro.NormErr('var', T, exFunH, 'Linf', i);
end% for

PrintTable.L2   = errL2;
PrintTable.Linf = errLinf;

PrintTable
end

