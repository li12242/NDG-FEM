function c = StyDiffSolver(mesh, c, q, TOLERR)
% 1D Steady diffusion solver
% Purpose: solver

resC = c; resQ = q;
err = 1;

while err>TOLERR
    [c, q] = StyDiffItera(mesh, resC, resQ);
    err = sum( (c(:) - resC(:)).^2 + (q(:) - resQ(:)).^2 );
    resC = c; 
    resQ = q;
    fprintf('error = %f\n', err)
end% while
end% func