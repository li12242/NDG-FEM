function outfile = CreateOutputFile(fileName, mesh)
% create netcdf file

% create nc file
ncid = netcdf.create(fileName,'CLOBBER');

% define nc file dimension
[dim_coor_id, dim_time_id] = DefineDim(ncid, mesh);
% define nc file variable
[x_id, u_id, time_id] = DefVar(ncid, dim_coor_id, dim_time_id);


netcdf.endDef(ncid); netcdf.close(ncid);

% reorder output variable
outfile = setOutputVar(ncid, [dim_coor_id, dim_time_id],...
    {'coordinate', 'time'},...
    [x_id, u_id, time_id], ...
    {'x', 'vairable', 'time'});
end% func

function [x_id, u_id, time_id] = DefVar(ncid, dim_coor_id, dim_time_id)
% define variable type & size
% Input:
%   ncid    - handle of netcdf file
%   dim_coor_id - indictor of coordinate dimension
%   dim_time_id - indictor of coordinate dimension
% Output:
%   x_id    - variable indictor of node position
%   u_id    - variable indictor of scalar
%   time_id - variable indictor of time
% 
x_id = netcdf.defVar(ncid,'x','double',[dim_coor_id, dim_time_id]);
time_id = netcdf.defVar(ncid,'time','double', dim_time_id);
u_id = netcdf.defVar(ncid,'u','double',[dim_coor_id, dim_time_id]);

end% func

function [coor_id, time_id] = DefineDim(ncid, mesh)
% define dimension size
% Input:
%   ncid    - handle of netcdf file
%   mesh    - mesh object of line type
% Output:
%   coor_id - indictor of coordinate dimension
%   time_id - indictor of time dimension
coor_id = netcdf.defDim(ncid,'coor', numel(mesh.x));
time_id = netcdf.defDim(ncid,'time', netcdf.getConstant('NC_UNLIMITED'));
end%func

function outfile = setOutputVar(ncid, dimids, dimname, varids, varname)
% store netcdf file indictor
outfile.ncid = ncid;
outfile.dimid = dimids;
outfile.dimname = dimname;
outfile.varid = varids;
outfile.varName = varname;
end% func
