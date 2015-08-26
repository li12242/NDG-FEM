function testing_jacobiP(n)
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
% [~, p] = calllib('libpolylib', 'jacobiP', np, r, p, n, alpha, beta);
% unloadlibrary libpolylib
p = Polylib.JacobiP(r,alpha,beta,n);

%% NDG
p1 = JacobiP(r, alpha, beta, double(n));

%% compare

p1-p'
