function [ Noutput ] = accessOutputStepNumber( obj )
m = 1;
filename = [obj.getOption('outputNetcdfCaseName'), '.', ...
    num2str(m), '-', num2str(obj.Nmesh), '.nc'];

ncid = netcdf.open(filename);
timeDimId = netcdf.inqDimID( ncid, 'Nt' );
[ ~, Noutput ] = netcdf.inqDim(ncid, timeDimId);
netcdf.close(ncif);
end