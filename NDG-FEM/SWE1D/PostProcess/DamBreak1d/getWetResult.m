!ifort -c main.f90 mod_RiemannSolver.f90
!ifort main.o mod_RiemannSolver.o -o main

% input parameters
chlen = 1000;
damposition = 500;
gravity = 9.81;
mCells = 500;
Tol = 1e-6;
iterNum = 100;
finalTime = [4, 8, 12];
dl = 10; ul = 0;
dr = 2; ur = 0;


for itime = 1:numel(finalTime)
    % creat input files
    filename = 'datain.ini';
    fig = fopen(filename, 'w');
    fprintf(fig, '%f\n', chlen);
    fprintf(fig, '%f\n', damposition);
    fprintf(fig, '%f\n', gravity);
    fprintf(fig, '%d\n', mCells);
    fprintf(fig, '%f\n', Tol);
    fprintf(fig, '%d\n', iterNum);
    fprintf(fig, '%f\n', finalTime(itime));
    fprintf(fig, '%f\n', dl);
    fprintf(fig, '%f\n', ul);
    fprintf(fig, '%14.12f\n', dr);
    fprintf(fig, '%f\n', ur);
    fclose(fig);
    
    !./main datain.ini
    data = load('result.out');
    x1 = data(:, 1); h1 = data(:, 2); u1 = data(:, 3); q1 = data(:, 4);
    save(['DamBreakWet', num2str(finalTime(itime)), 'mat'], 'x1', 'h1', 'u1', 'q1');
end
