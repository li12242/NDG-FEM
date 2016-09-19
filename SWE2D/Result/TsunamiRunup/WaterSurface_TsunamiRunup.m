function WaterSurface_TsunamiRunup
filename = 'output_ch5-7-9.txt';
fp = fopen(filename);
fgetl(fp);
ndata = fscanf(fp, '%f %f %f %f', [4, inf]);
fclose(fp);

xp = [4.521, 4.521, 4.521];
yp = [1.196, 1.696, 2.196];

% result
meshtype = 'quad';
filename = {'SWE2D.nc'};
order    = 1;
Postpro  = Utilities.PostProcess.Postprocess(filename, meshtype, order);
fileID   = 1;
time     = Postpro.NcFile(fileID).GetVarData('time');
bot      = Postpro.NcFile(fileID).GetVarData('bot');
nt       = numel(time);

% plot water surface at spicific position
for i = 1:3
    subplot(3,1,i); hold on;
    plot(ndata(1,:),ndata(1+i,:)/100,'r-');
    hp = zeros(nt, 1);
    for t = 1:nt
        h     = Postpro.NcFile(fileID).GetTimeVarData('h', t);
        eta   = h+bot;
        hp(t) = Postpro.Interp2D(eta, xp(i), yp(i), fileID);
        fprintf('Processing...%f\n', t/nt/3+(i-1)/3);
    end% for
    plot(time-0.65, hp, 'b.-');
end% func
end% func