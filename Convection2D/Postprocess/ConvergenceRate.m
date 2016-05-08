function ConvergenceRate

degree = 1:6;
ele = [40, 60, 80, 100];

parfor ideg = 1:numel(degree)
    e2 = zeros(numel(ele), 1);
    eInf = zeros(numel(ele), 1);
    
    for ine = 1: numel(ele)
        [mesh, var] = Convection2DSetUp(degree(ideg), ele(ine));
        
        [ev] = exactSolution(mesh.x, mesh.y);
        
        e2(ine) = L2(var, ev);
        eInf(ine) = Linf(var, ev);
        
    end% func
    
    rate2 = calConvRate(e2, ele);
    rateInf = calConvRate(eInf, ele);
    
    fig = fopen(['n', num2str(degree(ideg)), '.txt'],'w');
    fprintf('========================= n = %d =========================\n',degree(ideg));
    fprintf('----------------------- rate of L2 -----------------------\n');
    printResult(fig, rate2, e2, ele);
    fprintf('---------------------- rate of Linf ----------------------\n');
    printResult(fig, rateInf, eInf, ele);
    fclose(fig);
    
    loglog(1./ele, e2, 'ro-'); hold on;
    loglog(1./ele, eInf, 'bo-');
    t = legend('$L_2$', '$L_{\infty}$');
    set(t, 'box', 'off', 'Interpreter', 'Latex');
    xlabel('$\Delta x$', 'Interpreter', 'Latex');
    ylabel('Error', 'Interpreter', 'Latex');

end% func

end% func

function printResult(fig, rate, error, nele)
ne = numel(nele);

fprintf(fig, '|nele \t| L \t\t| Rate |\n');
formatStr = '|%d \t|%8.2e \t|%4.2f \t|\n';
        
for ie = 1:ne
    fprintf(fig, formatStr, nele(ie), error(ie), rate(ie));
end

p = polyfit(log(1./nele'), log(error), 1);
fprintf(fig, '|Fitted, \t|\\ \t\t\t|%4.2f \t|\n', p(1));
end% func

function rate = calConvRate(err, nele)
ne = numel(nele);
dx = 1./nele;
if ~iscolumn(dx)
    dx = dx';
end

dx = repmat(dx, 1, size(err, 2));
rate = zeros( size(err) );

for i = 2:ne
    rate(i, :) = log2( err(i-1,:)./err(i,:) )./log2( dx(i-1,:)./dx(i,:) );
end
end

function [ev] = exactSolution(x, y)

sigma = 125*1e3/33^2; 
xc = 0; yc = 3/5;
ev = exp(-sigma.*( (x - xc).^2 + (y - yc).^2) );

end% func

function err = L2(y, ey)
err = sqrt( sum2((y - ey).^2)./numel(y) );
end% func

function err = Linf(y, ey)
err = max2( abs(y - ey) );
end

function sumM = sum2(M)
sumM = sum(sum(M));
end

function maxM = max2(M)
maxM = sum(sum(M));
end% func