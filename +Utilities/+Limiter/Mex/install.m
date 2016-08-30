%% Limiter for 1-d
outdir = '../+Limiter1d';
src{1} = 'Minmod1d_Mex.c';
src{2} = 'Limiter.c';

mex('-O', src{:}, '-outdir', outdir);

%% Limiter for 2-d
outdir = '../+Limiter2d';
src    = 'SL2d_Mex.c';
mex('-O', src, '-outdir', outdir);

src    = 'SLLoc2d_Mex.c';
mex('-O', src, '-outdir', outdir);

src{1}    = 'TVB2d_Mex.c';
src{2}    = 'MatrixUtilities.c';
src{3}    = 'BasicFunction.c';
mex('-O', src{:}, '-outdir', outdir);

src{1}    = 'HWENO2d_Mex.c';
src{2}    = 'MatrixUtilities.c';
src{3}    = 'BasicFunction.c';
mex('-O', src{:}, '-outdir', outdir);

src{1}    = 'VB2d_Mex.c';
src{2}    = 'MatrixUtilities.c';
src{3}    = 'BasicFunction.c';
mex('-O', src{:}, '-outdir', outdir);
