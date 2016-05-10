function StoreVar(ncfile, var, time, contour)

% ncid = ncfile.ncid;
ncid = netcdf.open(ncfile.filename,'WRITE');
% netcdf.putVar(ncid,varid,data)
nx = numel(var);
% time
time_id = ncfile.varid(4);
netcdf.putVar(ncid, time_id, contour, 1, time)
% var
var_id = ncfile.varid(3);
netcdf.putVar(ncid, var_id, [0, contour], [nx, 1], var(:))

netcdf.close(ncid);
end% func