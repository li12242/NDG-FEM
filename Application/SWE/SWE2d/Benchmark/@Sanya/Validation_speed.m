function Validation_speed( obj, id )

MeasureSpeed =[ fileparts( mfilename('fullpath') ), '/tide/MeasureSpeed.txt'];

Mspeed = load(MeasureSpeed);

time = ncread('Sanya2k_0613.1-1.nc', 'time');
Cspeed = obj.PostProcess_speed;

x1 = Mspeed(:,1);
y1 = Mspeed(:,id+1);
x2 = time(167:316,1);
y2 = Cspeed(id,167:316)';

plot(x1,y1,'or',x2,y2,'-b');
axis([100000 200000 0 1.0]);

end

