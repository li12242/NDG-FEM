filename = 'SWE1D_LakeAtRest.nc';
filename1 = 'SWE1D_LakeAtRest1.nc'; % h-type refinement
time = ncread(filename, 'time');
time1 = ncread(filename1, 'time');
% timestep = numel(time)
x = ncread(filename, 'x');
startStep = 1; itime = startStep;

% bottom topography
bedElevation = zeros(size(x));
a = 1.2; rm = 0.4; r = abs(x - 0.5);
index = (r < rm);
bedElevation(index) = a*exp(-0.5./(rm.^2 - r(index).^2))./exp(-0.5./rm^2);

% exact solution of water depth
etaExact = ones(size(x));
Dryindex = etaExact < bedElevation;
etaExact(Dryindex) = bedElevation(Dryindex);

errQ = zeros(size(time)); errH = zeros(size(time));
errQ1 = zeros(size(time1)); errH1 = zeros(size(time1));

for itime = 1:numel(time)
    q = ncread(filename, 'q', [1, itime],[inf, 1]);
    h = ncread(filename, 'h', [1, itime],[inf, 1]);
    eta = h + bedElevation;
    
    errH(itime) = max(abs(eta - etaExact));
    errQ(itime) = max(abs(q));
    
    fprintf('Processing: %f ...\n', itime/numel(time))
end

for itime = 1:numel(time1)
    h1 = ncread(filename1, 'h', [1, itime],[inf, 1]);
    q1 = ncread(filename1, 'q', [1, itime],[inf, 1]);
    eta1 = h1 + bedElevation;

    errH1(itime) = max(abs(eta1 - etaExact));
    errQ1(itime) = max(abs(q1));
    
    fprintf('Processing: %f ...\n', itime/numel(time))
end

figure;
plot(time, errH, 'b.-'); hold on;
plot(time1, errH1, 'r.-')
xlabel('t'); ylabel('$\|\eta - \eta_0 \|_{\infty}$')
t = legend('convectional wet/dry treatment','h-type wet/dry treatment');
set(t, 'Box', 'off');

figure;
plot(time, errQ, 'b.-'); hold on;
plot(time1, errQ1, 'r.-')
xlabel('t'); ylabel('$\| q \|_{\infty}$')
t = legend('convectional wet/dry treatment','h-type wet/dry treatment');
set(t, 'Box', 'off');