function [e2, einf, r2, rinf] = NumAccuracy
%% Parameters
meshtype = 'tri';
filename = ['SWE2D_',meshtype,'_'];
degree   = 1;
ele      = [80, 100, 120, 160, 200];

% period
g     = 9.81;
alpha = 1.6*1e-7;
w     = sqrt(8*g*alpha);
T     = 2*pi/w;
time  = T; % spicific time

%% Loop to calculate error and convergence rate
e2    = zeros(numel(ele), 3); % error in L2 norm of h, qx, qy
einf  = zeros(numel(ele), 3);
vname = {'h', 'qx', 'qy'};
for ideg = 1:numel(degree)
    for ie = 1:numel(ele)
        % get the result at spicific time
        casename    = [filename, num2str(ele(ie)), '.nc'];
        [x, y, h, qx, qy] = GetResult(casename, time);
        % get exact result at spicific time
        [he, qxe, qye] = ParabolicBowlExtSol(x, y, time);
        % calculate the error in L1 and L2 norm
        [e2(ie, :), einf(ie, :)] = CalNormErr(h, qx, qy, he, qxe, qye);
    end% for
    % calculate the convergence rate
    r2   = CalConvReate(ele, e2);
    rinf = CalConvReate(ele, einf);
    
    for iv = 1:3
        fig = fopen([vname{iv},'.txt'],'w');
        fprintf('====================== %s ============================\n', vname{iv});
        PrintTable(fig, r2(:,iv), e2(:,iv), rinf(:,iv), einf(:,iv),...
            ele, degree(ideg), meshtype);
        fclose(fig);
    end
end% for
end% func

%% Print result
function PrintTable(fig, rate2, e2, rateInf, einf, nele, deg, type)
ne = numel(nele);

switch type
    case 'tri'
        dofs = nele.^2*2 * (deg+1)*(deg+2)/2;
    case 'quad'
        dofs = nele.^2 * (deg+1)*(deg+1);
end% switch

fprintf(fig, '|nele \t|DOFs \t| L2 \t\t| Rate \t| Linf \t\t| Rate \t|\n');
fprintf(fig, '| --- 	| --- 	| --- 		| ---  | --- 		| --- 	|\n');
fprintf('nele \t DOFs \t L2 \t\t Rate \t Linf \t\t Rate \n');
fileFormatStr = '|%d \t|%d \t|%8.2e \t|%4.2f \t|%8.2e \t|%4.2f \t|\n';
formatStr = '%d \t %d \t %8.2e \t %4.2f \t %8.2e \t %4.2f\n';
for ie = 1:ne
    fprintf(fig, fileFormatStr, ...
        nele(ie), dofs(ie), e2(ie), rate2(ie), einf(ie), rateInf(ie));
    fprintf(formatStr, ...
        nele(ie), dofs(ie), e2(ie), rate2(ie), einf(ie), rateInf(ie));
end

p = polyfit(log(1./nele'), log(e2), 1);
q = polyfit(log(1./nele'), log(einf), 1);
fprintf(fig, '|Fitted \t|\\ \t|\\ \t|%4.2f \t|\\ \t|%4.2f \t|\n', p(1), q(1));
fprintf('Fitted \t \\ \t \\ \t %4.2f \t \\ \t %4.2f \t\n', p(1), q(1));
end% func

%% Convergence rate calculation
function rate = CalConvReate(ele, err)
ne = numel(ele);
dx = 1./ele;
if ~iscolumn(dx)
    dx = dx';
end% if

rate = zeros( size(err) );

for i = 2:ne
    rate(i, :) = log2( err(i-1,:)./err(i,:) )./log2( dx(i-1,:)./dx(i,:) );
end% for
end% func

%% Error calculation
function [e2, einf] = CalNormErr(h, qx, qy, he, qxe, qye)
e2      = zeros(1, 3); 
einf    = zeros(1, 3); 

e2(1)   = L2(h,  he);
e2(2)   = L2(qx, qxe);
e2(3)   = L2(qy, qye);
einf(1) = Linf(h,  he);
einf(2) = Linf(qx, qxe);
einf(3) = Linf(qy, qye);
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

%% Read result
function [x, y, h, qx, qy] = GetResult(casename, stime)
x        = ncread(casename, 'x');
y        = ncread(casename, 'y');
time     = ncread(casename, 'time');
ist      = numel( time );
terr     = abs(time(ist) - stime);
fprintf('Time Deviation: %f\n', terr);
[np, ne] = size(x);

% get result
h        = ncread(casename, 'h',  [1,1,ist], [np, ne, 1]);
qx       = ncread(casename, 'qx', [1,1,ist], [np, ne, 1]);
qy       = ncread(casename, 'qy', [1,1,ist], [np, ne, 1]);
end% func