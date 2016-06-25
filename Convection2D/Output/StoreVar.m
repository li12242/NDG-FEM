function StoreVar(ncfile, mesh, var, disFlag, time, contour)

% ncid = ncfile.ncid;
ncid = netcdf.open(ncfile.filename,'WRITE');
% netcdf.putVar(ncid,varid,data)
nx = mesh.nNode;
% time
time_id = ncfile.varid(5);
netcdf.putVar(ncid, time_id, contour, 1, time)
% var
var_id = ncfile.varid(3);
netcdf.putVar(ncid, var_id, [0, contour], [nx, 1], var(:))
% discontinuity detector
temp_id = ncfile.varid(4);
netcdf.putVar(ncid, temp_id, [0, contour], [mesh.nElement, 1], disFlag(:))

netcdf.close(ncid);
end% func