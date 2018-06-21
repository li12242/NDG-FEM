function initFromMesh( obj, mesh )

% set vtk output
if (mesh.type == enumMeshDim.One)
    % obj.vtkOutput = vtkOutput2d();
elseif (mesh.type == enumMeshDim.Two)
    obj.vtkOutput = VtkOutput2d(obj.casename, obj.Nfield, obj.timeInterval);
elseif (mesh.type == enumMeshDim.Three)
    obj.vtkOutput = VtkOutput3d(obj.casename, obj.Nfield, obj.timeInterval);
end

obj.vtkOutput.initFromMesh( mesh );

% define dimension
dimTime = NdgNcDim('Nt', 0);
dimK = NdgNcDim('K', mesh.K);
dimNp = NdgNcDim('Np', mesh.cell.Np);
dimNfield = NdgNcDim('Nvar', obj.Nfield);

% define variable
varTime = NdgNcVar('time', dimTime, enumNcData.NC_DOUBLE );
varField = NdgNcVar('fphys', [dimNp, dimK, dimNfield, dimTime], enumNcData.NC_DOUBLE);

% define file
obj.filename = [ obj.casename, '/', obj.casename, '.nc' ];
obj.ncfile = NdgNcFile( obj.filename, ...
    [dimTime, dimK, dimNp, dimNfield], [varTime, varField]);

% init file
obj.ncfile.defineIntoNetcdfFile();

% set properties
obj.timeVarableId = varTime.id;
obj.fieldVarableId = varField.id;
obj.mesh = mesh;

end