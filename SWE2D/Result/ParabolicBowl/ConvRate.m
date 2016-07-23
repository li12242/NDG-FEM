function ConvRate
meshtype = 'tri';
elenum   = [80, 100, 120, 160, 200];
filename = cell(numel(elenum), 1);
for i =1:numel(elenum)
    filename{i} = ['SWE2D_', meshtype, '_', num2str(elenum(i)), '.nc'];
end

PostproSWE2D = Utilities.PostProcess.Postprocess(filename, meshtype, 1);

% time period
g     = 9.81;
alpha = 1.6*1e-7;
w     = sqrt(8*g*alpha);
T     = 2*pi/w;

errL2        = zeros(PostproSWE2D.nfiles, 1);
errLinf      = zeros(PostproSWE2D.nfiles, 1);

% error for H
PrintTable   = table;
PrintTable.M = elenum';
PrintTable.dofs = PostproSWE2D.GetDofs;
for i =1:PostproSWE2D.nfiles
    errL2(i)   = PostproSWE2D.NormErr('h', T, @ParabolicBowlExtDepth, 'L2', i);
    errLinf(i) = PostproSWE2D.NormErr('h', T, @ParabolicBowlExtDepth, 'Linf', i);
end% for
% convergence rate for H
a2   = PostproSWE2D.ConvRate('h', T, @ParabolicBowlExtDepth, 'L2');
ainf = PostproSWE2D.ConvRate('h', T, @ParabolicBowlExtDepth, 'Linf');

PrintTable.L2   = errL2;
PrintTable.a2   = a2;
PrintTable.Linf = errLinf;
PrintTable.ainf = ainf;

fprintf('The error for water depth');
PrintTable

% error for qx
for i =1:PostproSWE2D.nfiles
    errL2(i)   = PostproSWE2D.NormErr('qx', T, @ParabolicBowlExtQx, 'L2', i);
    errLinf(i) = PostproSWE2D.NormErr('qx', T, @ParabolicBowlExtQx, 'Linf', i);
end% for
% convergence rate for qx
a2   = PostproSWE2D.ConvRate('qx', T, @ParabolicBowlExtQx, 'L2');
ainf = PostproSWE2D.ConvRate('qx', T, @ParabolicBowlExtQx, 'Linf');

PrintTable.L2   = errL2;
PrintTable.a2   = a2;
PrintTable.Linf = errLinf;
PrintTable.ainf = ainf;

fprintf('The error for flow flux x');
PrintTable

% error for qy
for i =1:PostproSWE2D.nfiles
    errL2(i)   = PostproSWE2D.NormErr('qy', T, @ParabolicBowlExtQy, 'L2', i);
    errLinf(i) = PostproSWE2D.NormErr('qy', T, @ParabolicBowlExtQy, 'Linf', i);
end% for
% convergence rate for qy
a2   = PostproSWE2D.ConvRate('qy', T, @ParabolicBowlExtQy, 'L2');
ainf = PostproSWE2D.ConvRate('qy', T, @ParabolicBowlExtQy, 'Linf');

PrintTable.L2   = errL2;
PrintTable.a2   = a2;
PrintTable.Linf = errLinf;
PrintTable.ainf = ainf;

fprintf('The error for flow flux y');
PrintTable
end% func