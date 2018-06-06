%> @brief Get the total number of output step in the output NetCDF file
function [ Noutput ] = accessOutputStepNumber( obj )

outputStepNum = zeros( obj.Nmesh, 1 );
for m = 1:obj.Nmesh
    ncid = netcdf.open( obj.outputFile{m} );
    timeDimId = netcdf.inqDimID( ncid, 'Nt' );
    [ ~, outputStepNum(m) ] = netcdf.inqDim(ncid, timeDimId);
    netcdf.close(ncid);
end

% check whether the output step numbers in each file are equal
if max( outputStepNum ) ~= min( outputStepNum )
    msgID = [mfilename, ':TheOutputStepIsNotEqual'];
    msgtext = 'The output step in all files are not equal.';
    throw( MException(msgID, msgtext) );
end

Noutput = outputStepNum(1);
end