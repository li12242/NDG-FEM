function ParabolicBowlCovRate
order = [1];
nele = [50, 100, 200, 400, 600, 800, 1000, 1500, 2000];

% coefficients
a = 3000; h0 = 10; g = 9.81; B = 5; w = sqrt(2*g*h0)./a;
T = 2*pi*a/sqrt(2*g*h0); FinalTime = T*0.5;

errorType = 2; % L2 or L2 & Linf

for ideg = 1:numel(order)
    errorH = zeros(numel(nele) ,errorType);
    errorQ = zeros(numel(nele) ,errorType);
    
    for ine = 1:numel(nele)
        filename = ['SWE1D_', num2str(order(ideg)), '_', num2str(nele(ine)), '.nc'];

        % read results
        x = ncread(filename, 'x'); time = ncread(filename, 'time');
        % bottom topography
        bed = h0.*(x.^2./a^2 - 1);
        % get output step & results
        [~, index] = min( abs(time - FinalTime) );
        h = ncread(filename, 'h', [1, index],[inf, 1]);
        q = ncread(filename, 'q', [1, index],[inf, 1]);
        
        % get exact results
        [he, qe, ~] = exactParabolicBowlSolution(FinalTime, x, bed);
        
%         plot(x, h+bed, '-b.', x, he+bed, 'r-.', x, bed, 'k');
%         plot(x, h-he, '-b.');
        
        % get errors
        switch errorType
            case 1
                errorH(ine, 1) = L2(h, he);
                errorQ(ine, 1) = L2(q, qe);
                
            case 2
                errorH(ine, 1) = L2(h, he);
                errorH(ine, 2) = Linf(h, he);
                errorQ(ine, 1) = L2(h, he);
                errorQ(ine, 2) = Linf(h, he);
        end
        
    end% for
    
    rateH = calConvRate(errorH, nele);
    rateQ = calConvRate(errorQ, nele);
    
    fprintf('========================= n = %d =========================\n',order(ideg));
    fprintf('----------------------- rate of h -----------------------\n');
    printResult(rateH, errorH, nele, errorType);
    fprintf('----------------------- rate of q -----------------------\n');
    printResult(rateQ, errorQ, nele, errorType);
    fprintf('\n');
end% for

end% func

function printResult(rate, error, nele, errorType)
ne = numel(nele); figure
switch errorType
    case 1
        fprintf('nele, \t L2\n');
        formatStr = '%d, \t %f\n';
        loglog(1./nele, error, 'o-');
        ylabel('$L_2$', 'Interpreter', 'Latex');
    case 2
        fprintf('nele, \t L2, \t Linf\n')
        formatStr = '%d, \t %f, \t %f\n';
        subplot(2, 1, 1); loglog(1./nele, error(:, 1), 'o-');
        ylabel('$L_2$', 'Interpreter', 'Latex');
        subplot(2, 1, 2); loglog(1./nele, error(:, 2), 'o-');
        ylabel('$L_{\infty}$', 'Interpreter', 'Latex');
end
for ie = 1:ne
    fprintf(formatStr, nele(ie), rate(ie, :));
end
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

function err = L2(y, ey)
err = sqrt( sum((y - ey).^2)./numel(y) );
end% func

function err = Linf(y, ey)
err = max( abs(y - ey) );
end