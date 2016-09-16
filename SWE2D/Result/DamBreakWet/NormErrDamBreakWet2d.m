function NormErrDamBreakWet2d
%NORMERRDAMBREAKDRY2D Summary of this function goes here
%   Detailed explanation goes here

%% parameters
meshtype = 'quad';
N        = 1;
filename = {'SWE2D.nc'};
Postpro = Utilities.PostProcess.Postprocess(filename, meshtype, N);

T = 20;
exFunH  = @DamBreakWet2d_H;
% exFunQx = @DamBreakWet2d_Qx;

%% table for norm error
errL2        = zeros(Postpro.nfiles, 1);
errLinf      = zeros(Postpro.nfiles, 1);

PrintTable   = table;
% PrintTable.M = elenum';
PrintTable.dofs = Postpro.GetDofs;

for i =1:Postpro.nfiles
    errL2(i)   = Postpro.NormErr('h', T, exFunH, 'L2', i);
    errLinf(i) = Postpro.NormErr('h', T, exFunH, 'Linf', i);
end% for

PrintTable.H_L2   = errL2;
PrintTable.H_Linf = errLinf;

% for i =1:Postpro.nfiles
%     errL2(i)   = Postpro.NormErr('qx', T, exFunQx, 'L2', i);
%     errLinf(i) = Postpro.NormErr('qx', T, exFunQx, 'Linf', i);
% end% for
% 
% PrintTable.Qx_L2   = errL2;
% PrintTable.Qx_Linf = errLinf;

PrintTable
end

