function ConvergenceRate

degree = 1:8;
ele = [40, 60, 80];

for ideg = 1:numel(degree)
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
    
    fprintf('========================= n = %d =========================\n',order(ideg));
    fprintf('----------------------- rate of L2 -----------------------\n');
    printResult(rate2, e2, ele, 'L_2');
    fprintf('---------------------- rate of LInf ----------------------\n');
    printResult(rateInf, eInf, ele, 'L_{\infty}');
    fprintf('\n');
    
    loglog(1./nele, e2, 'ro-'); hold on;
    loglog(1./nele, eInf, 'bo-');
    t = legend('L_2', 'L_{\infty}');
    set(t, 'box', 'off', 'Interpreter', 'Latex');
    ylabel('Error', 'Interpreter', 'Latex');

end% func

end% func

function printResult(rate, error, nele, errName)
ne = numel(nele); figure

fprintf('|nele \t| L \t\t| Rate |\n');
formatStr = '|%d \t|%8.2e \t|%4.2f \t|\n';
        
for ie = 1:ne
    fprintf(formatStr, nele(ie), error(ie), rate(ie));
end

p = polyfit(log(1./nele'), log(error), 1);
fprintf('|Fitted, \t|\\ \t\t\t|%4.2f \t|\n', p(1));
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