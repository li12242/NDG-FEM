function initFromMesh( obj, mesh )

% define dimension
dimTime = NdgNcDim('Nt', 0);
dimK = NdgNcDim('K', mesh.K);
dimNp = NdgNcDim('Np', mesh.cell.Np);
dimNfield = NdgNcDim('Nvar', obj.Nfield);

% define variable
varTime = NdgNcVar('time', dimTime, enumNcData.NC_DOUBLE );
varField = NdgNcVar('fphys', [dimNp, dimK, dimNfield, dimTime], enumNcData.NC_DOUBLE);

% define file
obj.filename = [ obj.casename, '.nc' ];
obj.ncfile = NdgNcFile( obj.filename, ...
    [dimTime, dimK, dimNp, dimNfield], [varTime, varField]);

% init file
obj.ncfile.defineIntoNetcdfFile();

% set properties
obj.timeVarableId = varTime.id;
obj.fieldVarableId = varField.id;
obj.mesh = mesh;
end