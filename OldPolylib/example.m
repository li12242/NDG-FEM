loadlibrary('libpolylib.dyld');
list = libfunctions('libpolylib','-full')
np = int32(5);
z = zeros(1,np); w = zeros(1,np);

[z, w] = calllib('libpolylib', 'zwgj', z, w,np,alpha,beta);

