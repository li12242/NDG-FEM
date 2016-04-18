function outfile = CreateOutputFile(mesh)

ncid = netcdf.create('SWE1D.nc','CLOBBER');

[coor_id, time_id, boundary_id] = DefineDim(ncid, mesh);
[x_id, bx_id, h_id, q_id, time_id, wavespeed_id, status_id]...
    = DefVar(ncid, coor_id, time_id, boundary_id);

netcdf.endDef(ncid);
netcdf.close(ncid);

outfile = setOutputVar(ncid, [coor_id, time_id, boundary_id],...
    {'coordinate', 'time', 'boundary'},...
    [x_id, h_id, q_id, time_id, wavespeed_id, status_id, bx_id], ...
    {'corrdinate', 'height', 'flux', 'time', 'wave speed', 'wet-dry status'});
end% func

function [x_id, bx_id, h_id, q_id, time_id, wavespeed_id, status_id] ...
    = DefVar(ncid, coor_id, time_id, boundary_id)
x_id = netcdf.defVar(ncid,'x','double',coor_id);
bx_id = netcdf.defVar(ncid, 'bx', 'double', boundary_id);
time_id = netcdf.defVar(ncid,'time','double', time_id);
h_id = netcdf.defVar(ncid,'h','double',[coor_id, time_id]);
q_id = netcdf.defVar(ncid,'q','double',[coor_id, time_id]);
wavespeed_id = netcdf.defVar(ncid,'wavespeed','double',time_id);
status_id = netcdf.defVar(ncid, 'status', 'double', [boundary_id, time_id]);

end% func

function [coor_id, time_id, boundary_id] = DefineDim(ncid,mesh)
coor_id = netcdf.defDim(ncid,'coor', numel(mesh.x));
boundary_id = netcdf.defDim(ncid, 'bound', numel(mesh.vmapM));
time_id = netcdf.defDim(ncid,'time', netcdf.getConstant('NC_UNLIMITED'));
end%func

function outfile = setOutputVar(ncid, dimids, dimname, varids, varname)
outfile.ncid = ncid;
outfile.dimid = dimids;
outfile.dimname = dimname;
outfile.varid = varids;
outfile.varName = varname;
end% func
