function genGoldenPolylib()
%> @brief Generate golden reference data for the Polylib gateway functions.
%>
%> Run this script BEFORE deleting the mex binaries (lib/*.mexw64): it
%> captures the mex outputs as the reference baseline for the pure-Matlab
%> port (see PolylibGatewayTest.m).
%>
%> NOTE: golden cases are limited to n<=12. For n>=13 the original C code
%> overflows an int32 factorial in the normalization factor (thirdParty/
%> Polylib/polylib.c), so the mex is systematically wrong there and its
%> output is not a valid reference. The Matlab port uses gamma() and is
%> exact for all n.
% ======================================================================
%> This function is part of the NDG-FEM software.
% ======================================================================

thisDir = fileparts(mfilename('fullpath'));
addpath( fullfile(thisDir, '..', '..', '..', 'lib') );

rng(0); % fixed seed for reproducibility
r = [ linspace(-1, 1, 17).'; -1 + 2*rand(25, 1) ];

alphaList = 0:12;
betaList = [0, 1];
nList = 0:12;
nCase = numel(alphaList) * numel(betaList) * numel(nList);

JP = repmat( struct('alpha', [], 'beta', [], 'n', [], 'r', [], 'P', []), nCase, 1 );
GP = repmat( struct('alpha', [], 'beta', [], 'n', [], 'r', [], 'dP', []), nCase, 1 );
c = 0;
for a = alphaList
    for b = betaList
        for n = nList
            c = c + 1;
            JP(c).alpha = a; JP(c).beta = b; JP(c).n = n;
            JP(c).r = r;  JP(c).P = JacobiP(r, a, b, n);
            GP(c).alpha = a; GP(c).beta = b; GP(c).n = n;
            GP(c).r = r;  GP(c).dP = GradJacobiP(r, a, b, n);
        end
    end
end

npList = 1:13;
Z = repmat( struct('np', [], 'z', [], 'w', []), numel(npList), 1 );
for k = npList
    [z, w] = zwglj(k);
    Z(k).np = k; Z(k).z = z; Z(k).w = w;
end

golden = struct( ...
    'jacobiP',     {JP}, ...
    'gradJacobiP', {GP}, ...
    'zwglj',       {Z}  );
golden.meta = struct( ...
    'createdBy', 'genGoldenPolylib (mex baseline)', ...
    'matlabVersion', version, ...
    'arch', computer('arch') );

matFile = fullfile(thisDir, 'polylib_golden.mat');
save(matFile, 'golden');
fprintf('Golden reference saved to %s\n', matFile);
fprintf('  jacobiP cases:     %d\n', numel(golden.jacobiP));
fprintf('  gradJacobiP cases: %d\n', numel(golden.gradJacobiP));
fprintf('  zwglj cases:       %d\n', numel(golden.zwglj));
end% func
