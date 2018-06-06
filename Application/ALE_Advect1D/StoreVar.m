function StoreVar(fileName, ncfile, x, u, time, contour)

% ncid = ncfile.ncid;
ncid = netcdf.open(fileName,'WRITE');

nx = numel(x);
% time
time_id = ncfile.varid(3);
netcdf.putVar(ncid, time_id, contour, 1, time)
% x
x_id = ncfile.varid(1);
netcdf.putVar(ncid, x_id, [0, contour], [nx, 1], x(:))
% u
u_id = ncfile.varid(2);
netcdf.putVar(ncid,u_id, [0, contour], [nx, 1],u(:))

netcdf.close(ncid);
end% func