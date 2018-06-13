function Validation_level( obj, id1, id2 )
MeasureLevel =[ fileparts( mfilename('fullpath') ), '/tide/MeasureLevel.txt'];
% MeasureLevel ='I:\MeasureLevel.txt';
Mlevel = load(MeasureLevel);

time = ncread('Sanya2k_0613.1-1.nc', 'time');
Clevel = obj.PostProcess_level;

x1 = Mlevel(:,1);
y1 = Mlevel(:,id1+1);
x2 = time(167:316,1);
y2 = Clevel(id2,167:316)';

plot(x1,y1,'or',x2,y2,'-b'); 
axis([100000 200000 -0.8 1.2]);
end