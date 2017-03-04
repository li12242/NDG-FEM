function NormErr_ObliqueHydraulicJump
%NORMERR_OBLIQUEHYDRAULICJUMP Summary of this function goes here
%   Detailed explanation goes here

meshtype = 'quad';
N        = 1;
filename = {'ObliqueHydraulicJump_VBHWENO_quad.nc', ...
    'ObliqueHydraulicJump_VBHWENO_tri'};
Postpro = Utilities.PostProcess.Postprocess(filename, meshtype, N);
fileID = 1;

time = Postpro.NcFile(fileID).GetVarData('time');
nt = numel(time);
exFunH  = @ObliqueHydraulicJump_H;

%% table for norm error
errL2        = zeros(nt, Postpro.nfiles);
errLinf      = zeros(nt, Postpro.nfiles);

for i =1:nt
    for ifile = 1:Postpro.nfiles
        errL2(i, ifile)   = Postpro.NormErr('h', time(i), exFunH, 'L2', fileID);
        errLinf(i, ifile) = Postpro.NormErr('h', time(i), exFunH, 'Linf', fileID);
        fprintf('Processing %f...\n', i/nt);
    end
end% for

for ifile = 1:Postpro.nfiles
    plot(time, errL2(:, ifile), 'rs-'); hold on;
    plot(time, errLinf(:, ifile), 'bo-');
end% for

legend({'L2', 'Linf'}, 'box', 'off')
end

