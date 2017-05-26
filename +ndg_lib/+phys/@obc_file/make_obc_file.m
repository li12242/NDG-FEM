function make_obc_file( obj, filename, time, vert, f_extQ )
%MAKE_OBC_FILE Summary of this function goes here
%   Detailed explanation goes here

[Nv, Nt, Nfield] = size(f_extQ);
time_d = ndg_utility.nc.nc_dim('Nt', Nt);
nfld_d = ndg_utility.nc.nc_dim('Nfield', Nfield);
nv_d = ndg_utility.nc.nc_dim('Nv', Nv);

time_v = ndg_utility.nc.nc_var('time', time_d, 'double');
vert_v = ndg_utility.nc.nc_var('vert', nv_d, 'int');
fext_v = ndg_utility.nc.nc_var('f_extQ', [nv_d, nfld_d, time_d], 'double');

obcfile = ndg_utility.nc.nc_file(filename, ...
    [nv_d, time_d, nfld_d], [time_v, vert_v, fext_v]);

obcfile.ncid = netcdf.create(filename,'CLOBBER');
for n = 1:numel(obcfile.dims)
    obcfile.dims(n).define_in_ncfile(obcfile.ncid);
end

for n = 1:numel(obcfile.vars)
    obcfile.vars(n).define_in_ncfile(obcfile.ncid);
end
netcdf.endDef(obcfile.ncid);

% put variables
netcdf.putVar(obcfile.ncid, obcfile.vars(1).id, time);
netcdf.putVar(obcfile.ncid, obcfile.vars(2).id, vert);
netcdf.putVar(obcfile.ncid, obcfile.vars(3).id, f_extQ);

netcdf.close(obcfile.ncid);
end

