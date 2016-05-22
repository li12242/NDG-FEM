function ConvergenceRate

degree = 1:2;
ele = [40, 60, 80];

for ideg = 1:numel(degree)
    e2 = zeros(numel(ele), 1);
    eInf = zeros(numel(ele), 1);
    
    for ine = 1: numel(ele)
        %[mesh, var] = Convection2DSetUp( degree(ideg), ele(ine));
        filename = ['Convection2D_', num2str(degree(ideg)),'_',num2str(ele(ine)),'.nc'];
        x = ncread(filename, 'x'); y = ncread(filename, 'y');
        time = ncread(filename, 'time');
        var = ncread(filename, 'var', [1, numel(time)],[inf, 1]);
        
        [ev] = exactSolution(x, y);
        
        e2(ine) = L2(var, ev);
        eInf(ine) = Linf(var, ev);
        
    end% func
    
    rate2 = calConvRate(e2, ele, degree(ideg));
    rateInf = calConvRate(eInf, ele, degree(ideg));
    
    fig = fopen(['n', num2str(degree(ideg)), '.txt'],'w');
    fprintf('========================= n = %d =========================\n',degree(ideg));
    printResult(fig, rate2, e2, rateInf, eInf, ele, degree(ideg));
    fclose(fig);
    
end% func

end% func

function printResult(fig, rate2, e2, rateInf, einf, nele, deg)
ne = numel(nele);

type = 'tri';
switch type
    case 'tri'
        dofs = nele.^2*2 * (deg+1)*(deg+2)/2;
    case 'quad'
        dofs = nele.^2 * (deg+1)*(deg+1);
end

fprintf(fig, '|nele \t|DOFs \t| L2 \t\t| Rate \t| Linf \t\t| Rate \t|\n');
fprintf(fig, '| --- 	| --- 	| --- 		| ---  | --- 		| --- 	|\n');
fprintf('nele \t DOFs \t L2 \t\t Rate \t Linf \t\t Rate \n');
fileFormatStr = '|%d \t|%d \t|%8.2e \t|%4.2f \t|%8.2e \t|%4.2f \t|\n';
formatStr = '%d \t %d \t %8.2e \t %4.2f \t %8.2e \t %4.2f\n';
for ie = 1:ne
    fprintf(fig, fileFormatStr, nele(ie), dofs(ie), e2(ie), rate2(ie), einf(ie), rateInf(ie));
    fprintf(formatStr, nele(ie), dofs(ie), e2(ie), rate2(ie), einf(ie), rateInf(ie));
end

p = polyfit(log(1./nele'), log(e2), 1);
q = polyfit(log(1./nele'), log(einf), 1);
fprintf(fig, '|Fitted \t|\\ \t|\\ \t|%4.2f \t|\\ \t|%4.2f \t|\n', p(1), q(1));
fprintf('Fitted \t \\ \t \\ \t %4.2f \t \\ \t %4.2f \t\n', p(1), q(1));
end% func

function rate = calConvRate(err, nele, deg)
ne = numel(nele);
dx = 1./nele;
if ~iscolumn(dx)
    dx = dx';
end

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
maxM = max(max(M));
end% func