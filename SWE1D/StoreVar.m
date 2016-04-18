function StoreVar(ncfile, h, q, time, lamda, contour)

% ncid = ncfile.ncid;
ncid = netcdf.open('SWE1D.nc','WRITE');
% netcdf.putVar(ncid,varid,data)
nx = numel(h);

% h
h_id = ncfile.varid(2);
netcdf.putVar(ncid,h_id, [0, contour], [nx, 1], h(:))
% q
q_id = ncfile.varid(3);
netcdf.putVar(ncid,q_id, [0, contour], [nx, 1],q(:))
% time
time_id = ncfile.varid(4);
netcdf.putVar(ncid, time_id, contour, 1, time)
% wavespeed
wavespeed_id = ncfile.varid(5);
netcdf.putVar(ncid,wavespeed_id, contour, 1, lamda)

% % wet/dry status
% nboundary = numel(status);
% status_id = ncfile.varid(6);
% netcdf.putVar(ncid, status_id, [0, contour], [nboundary, 1], status)

netcdf.close(ncid);
end% func