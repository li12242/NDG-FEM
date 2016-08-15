outdir = '../+Limiter1d';
src{1} = 'Minmod1d_Mex.c';
src{2} = 'Limiter.c';

mex('-O', src{:}, '-outdir', outdir);