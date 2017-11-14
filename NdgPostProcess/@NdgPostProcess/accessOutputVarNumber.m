function [ Nvar ] = accessOutputVarNumber( obj )
%ACCESSOUTPUTVARNUMBER Summary of this function goes here
%   Detailed explanation goes here

outputVarNum = zeros( obj.Nmesh, 1 );
for m = 1:obj.Nmesh
    ncid = netcdf.open( obj.outputFile{m} );
    timeDimId = netcdf.inqDimID( ncid, 'Nvar' );
    [ ~, outputVarNum(m) ] = netcdf.inqDim(ncid, timeDimId);
    netcdf.close(ncid);
end

% check whether the output step numbers in each file are equal
if max( outputVarNum ) ~= min( outputVarNum )
    msgID = [mfilename, ':TheOutputVariableNumberIsNotEqual'];
    msgtext = 'The output variable number in all files are not equal.';
    throw( MException(msgID, msgtext) );
end

Nvar = outputVarNum(1);
end

