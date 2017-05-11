function declare_file( obj )
%DECLARE_FILE Create a new NetCDF file and define contained variables
%   Detailed explanation goes here

obj.ncid = netcdf.create(obj.name,'CLOBBER');

for n = 1:numel(obj.dims)
    obj.dims(n).define_in_ncfile(obj.ncid);
end

for n = 1:numel(obj.vars)
    obj.vars(n).define_in_ncfile(obj.ncid);
end
netcdf.endDef(obj.ncid);
obj.isopen = true;

end

