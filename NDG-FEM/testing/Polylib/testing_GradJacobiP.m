function testing_GradJacobiP(n)
r = Polylib.zwglj(n+1);
alpha=1; beta=1;
%% poylib
% hpath = fullfile('/Users/mac/Documents/Model/NDG-FEM/+Polylib/libpolylib.h');
% libpath = fullfile('/Users/mac/Documents/Model/NDG-FEM/+Polylib/libpolylib.dyld');
% loadlibrary(libpath,hpath);
% 
% np = int32(numel(r));
% n = int32(n);
% % alpha = 0; beta =0;
% p = zeros(1,np);
% 
% 
% [~, p] = calllib('libpolylib', 'GradjacobiP', np, r, p, n, alpha, beta);
% unloadlibrary libpolylib
p = Polylib.GradJacobiP(r,alpha,beta,n);

%% NDG
p1 = GradJacobiP(r, alpha, beta, double(n));