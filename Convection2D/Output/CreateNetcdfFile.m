function outfile = CreateNetcdfFile(mesh, filename)

ncid = netcdf.create(filename,'CLOBBER');

[node_id, ele_id, time_id] = DefineDim(ncid, mesh);
[x_id, y_id, var_id, temp_id, time_id] = DefVar(ncid, node_id, ele_id, time_id);

netcdf.endDef(ncid);
netcdf.close(ncid);

outfile = setOutputVar(filename, [node_id, ele_id, time_id],...
    {'node', 'ele', 'time'},...
    [x_id, y_id, var_id, temp_id, time_id], ...
    {'x', 'y', 'scalar', 'temp', 'time'});

ncid = netcdf.open(filename,'WRITE');
netcdf.putVar(ncid, x_id, 0, numel(mesh.x), mesh.x);
netcdf.putVar(ncid, y_id, 0, numel(mesh.y), mesh.y);
netcdf.close(ncid);
end% func

function [x_id, y_id, var_id, temp_id, time_id] = DefVar(ncid, node_id, ele_id, time_id)
x_id = netcdf.defVar(ncid,'x','double', node_id);
y_id = netcdf.defVar(ncid,'y','double', node_id);
var_id = netcdf.defVar(ncid,'h','double',[node_id, time_id]);
temp_id = netcdf.defVar(ncid,'temp','double',[ele_id, time_id]);
time_id = netcdf.defVar(ncid,'time','double', time_id);
end% func

function [node_id, ele_id, time_id] = DefineDim(ncid,mesh)
node_id = netcdf.defDim(ncid,'node', mesh.nNode);
ele_id = netcdf.defDim(ncid, 'ele', mesh.nElement);
time_id = netcdf.defDim(ncid,'time', netcdf.getConstant('NC_UNLIMITED'));
end%func

function outfile = setOutputVar(filename, dimids, dimname, varids, varname)
outfile.filename = filename;
outfile.dimid = dimids;
outfile.dimname = dimname;
outfile.varid = varids;
outfile.varName = varname;
end% func