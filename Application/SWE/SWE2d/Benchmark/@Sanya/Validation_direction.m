function Validation_direction( obj, id )

MeasureDir =[ fileparts( mfilename('fullpath') ), '/tide/MeasureDirection.txt'];

Mdir = load(MeasureDir);

time = ncread('Sanya2k_0613.1-1.nc', 'time');
Cdir = obj.PostProcess_direction;

x1 = Mdir(:,1);
y1 = Mdir(:,id+1);
x2 = time(167:316,1);
y2 = Cdir(id,167:316)';

plot(x1,y1,'or',x2,y2,'-b');
axis([100000 200000 0 360]);


end

