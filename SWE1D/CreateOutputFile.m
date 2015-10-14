function outfile = CreateOutputFile(mesh)

ncid = netcdf.create('SWE1D.nc','CLOBBER');

[coor_id, time_id] = DefineDim(ncid, mesh);
[x_id, h_id, q_id, time_id, wavespeed_id] = DefVar(ncid, coor_id, time_id);

netcdf.endDef(ncid);
netcdf.close(ncid);

outfile = setOutputVar(ncid, [coor_id, time_id],...
    {'coordinate', 'time'},...
    [x_id, h_id, q_id, time_id, wavespeed_id], ...
    {'corrdinate', 'height', 'flux', 'time', 'wave speed'});
end% func

function [x_id, h_id, q_id, time_id, wavespeed_id] = DefVar(ncid, coor_id, time_id)
x_id = netcdf.defVar(ncid,'x','double',coor_id);
time_id = netcdf.defVar(ncid,'time','double', time_id);
h_id = netcdf.defVar(ncid,'h','double',[coor_id, time_id]);
q_id = netcdf.defVar(ncid,'q','double',[coor_id, time_id]);

wavespeed_id = netcdf.defVar(ncid,'wavespeed','double',time_id);

end% func

function [coor_id, time_id] = DefineDim(ncid,mesh)
coor_id = netcdf.defDim(ncid,'coor', numel(mesh.x));
time_id = netcdf.defDim(ncid,'time', netcdf.getConstant('NC_UNLIMITED'));
end%func

function outfile = setOutputVar(ncid, dimids, dimname, varids, varname)
outfile.ncid = ncid;
outfile.dimid = dimids;
outfile.dimname = dimname;
outfile.varid = varids;
outfile.varName = varname;
end% func
