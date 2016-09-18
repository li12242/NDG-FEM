function WaterSurface_TsunamiRunup
filename = 'output_ch5-7-9.xls';
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
x        = Postpro.NcFile(fileID).GetVarData('x');
y        = Postpro.NcFile(fileID).GetVarData('y');
bot      = Postpro.NcFile(fileID).GetVarData('bot');
interp   = TriScatteredInterp(x(:), y(:), bot(:), 'linear');

% plot water surface at spicific position
for i = 1:3
    subplot(3,1,i); hold on;
    plot(ndata(1,:),ndata(1+i,:),'r-');
    hp = zeros(size(time), 1);
    for t = 1:numel(time)
        h = Postpro.NcFile(fileID).GetTimeVarData('h', t);
        interp.V = h(:) + bot(:);
        hp(t) = interp(xp(i), yp(i));
    end% for
    plot(time, hp, 'b.-');
end% func
end% func