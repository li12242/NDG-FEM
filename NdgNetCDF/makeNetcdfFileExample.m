Nt = 20;
Nfield = 3;
Nv = 10;
time_d = NdgNcDim('Nt', Nt);
nfld_d = NdgNcDim('Nfield', Nfield);
nv_d = NdgNcDim('Nv', Nv);

time_v = NdgNcVar('time', time_d, NdgNcType.NC_DOUBLE);
vert_v = NdgNcVar('vert', nv_d, NdgNcType.NC_INT);
fext_v = NdgNcVar('f_extQ', [nv_d, nfld_d, time_d], NdgNcType.NC_DOUBLE);

filename = 'advUniformMesh.0-1.nc';
ncfile = NdgNcFile(filename, [nv_d, time_d, nfld_d], [time_v, vert_v, fext_v]);
ncfile.defineIntoNetcdfFile(); % define the NetCDF file;
ncfile.outputTotalVar( 1, 1:Nt );
