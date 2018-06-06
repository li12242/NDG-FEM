Nt = 20;
Nfield = 3;
Nv = 10;
time_d = NdgNcDim('Nt', Nt);
nfld_d = NdgNcDim('Nfield', Nfield);
nv_d = NdgNcDim('Nv', Nv);

time_v = NdgNcVar('time', time_d, NdgNcDataType.NC_DOUBLE);
vert_v = NdgNcVar('vert', nv_d, NdgNcDataType.NC_INT);
fext_v = NdgNcVar('f_extQ', [nv_d, nfld_d, time_d], NdgNcDataType.NC_DOUBLE);

filename = 'outputTest.0-1.nc';
ncfile = NdgNcOutputFile(filename, [nv_d, time_d, nfld_d], [time_v, vert_v, fext_v]);
ncfile.defineIntoNetcdfFile(); % define the NetCDF file;
ncfile.outputTotalVar( 1, 1:Nt );

outputfile = getOutputFile( NdgIOFileType.NetCDF, NdgIOIntervalType.DeltaTime, ...
    filename, [nv_d, time_d, nfld_d], [time_v, vert_v, fext_v], 0.1 );
