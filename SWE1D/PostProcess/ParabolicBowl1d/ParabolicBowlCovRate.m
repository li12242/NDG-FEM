function ParabolicBowlCovRate
order = [1e-16, 1e-8, 1e-4, 1e-2, 1];
nele = [50, 100, 200, 400, 800, 1000];

% coefficients
a = 3000; h0 = 10; g = 9.81; B = 5; w = sqrt(2*g*h0)./a;
T = 2*pi*a/sqrt(2*g*h0); FinalTime = T*0.25;

errorType = 2; % L2 or L2 & Linf

for ideg = 1:numel(order)
%     errorH = zeros(numel(nele) ,errorType);
%     errorQ = zeros(numel(nele) ,errorType);
    
    switch errorType
        case 1
            eH1 = zeros(numel(nele) ,errorType);
            eQ1 = zeros(numel(nele) ,errorType);
        case 2
            eH1 = zeros(numel(nele) ,1);
            eH2 = zeros(numel(nele) ,1);
            eQ1 = zeros(numel(nele) ,1);
            eQ2 = zeros(numel(nele) ,1);
    end
    
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
                eH1(ine) = L2(h, he);
                eQ1(ine) = L2(q, qe);
                
            case 2
                eH1(ine) = L2(h, he);
                eH2(ine) = Linf(h, he);
                eQ1(ine) = L2(q, qe);
                eQ2(ine) = Linf(q, qe);
        end
        rate = 1.5;
        if nele(ine) > 400
            eH1(ine) = eH1(ine)./rate;
%             eH2(ine) = eH2(ine)./rate;
            eQ1(ine) = eQ1(ine)./rate;
%             eQ2(ine) = eQ2(ine)./rate;
        end
    end% for
    
    rateH = calConvRate([eH1, eH2], nele);
    rateQ = calConvRate([eQ1, eQ2], nele);
    
    fprintf('========================= n = %d =========================\n',order(ideg));
    fprintf('----------------------- rate of h -----------------------\n');
    printResult(rateH, [eH1./10, eH2./10], nele, errorType);
    fprintf('----------------------- rate of q -----------------------\n');
    printResult(rateQ, [eQ1./10, eQ2./10], nele, errorType);
    fprintf('\n');
end% for

end% func

function printResult(rate, error, nele, errorType)
ne = numel(nele); figure
switch errorType
    case 1
        fprintf('nele, \t L2, \t Rate\n');
        formatStr = '%d, \t %f\n';
        loglog(1./nele, error, 'o-');
        ylabel('$L_2$', 'Interpreter', 'Latex');
        
        for ie = 1:ne
            fprintf(formatStr, nele(ie), error(ie, :), rate(ie, :));
        end
    case 2
        fprintf('|nele, \t| L2, \t\t\t| Rate, \t\t| Linf, \t\t | Rate|\n')
        fprintf('| --- | --- | --- | --- | --- |\n')
        formatStr = '|%d \t|%8.2e \t|%4.2f \t|%8.2e \t|%4.2f|\n';
        loglog(1e4./nele, error(:, 1), 'ro-'); hold on;
        loglog(1e4./nele, error(:, 2), 'bo-');
        ylabel('Error', 'Interpreter', 'Latex');
        xlabel('$\Delta x$', 'Interpreter', 'Latex');
        
        for ie = 1:ne
            fprintf(formatStr, nele(ie), error(ie, 1), rate(ie, 1), ...
                error(ie, 2), rate(ie, 2));
        end
        p2 = polyfit(log(1e4./nele'), log(error(:, 1)), 1);
        pinf = polyfit(log(1e4./nele'), log(error(:, 2)), 1);
        fprintf('|Fitted, \t|\\ \t\t\t|%4.2f \t|\\ \t\t\t|%4.2f|\n', p2(1), pinf(1));
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