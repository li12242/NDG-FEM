function matInitMatSolver( obj )

% setup the output file object
obj.outputNcFile = [];
for n = 1:obj.Nmesh
    [ ncfile ] = makeBasicOutputNetcdfFile( obj, n, obj.meshUnion(n) );
    obj.outputNcFile = [obj.outputNcFile, ncfile];
end

%> set the final time
obj.ftime = obj.getOption('finalTime');
%> choose the limiter 
obj.limiter = getLimiter( obj.getOption('LimiterType'), obj.meshUnion );
end% func

function [ ncfile ] = makeBasicOutputNetcdfFile( phys, ind, mesh )
dimTime = NdgNcDim('Nt', 0);
dimK = NdgNcDim('K', mesh.K);
dimNp = NdgNcDim('Np', mesh.cell.Np);
dimNfield = NdgNcDim('Nvar', phys.Nvar);

varTime = NdgNcVar('time', dimTime, NdgNcDataType.NC_DOUBLE );
varField = NdgNcVar('fphys', [dimNp, dimK, dimNfield, dimTime], NdgNcDataType.NC_DOUBLE);

filename = [ phys.getOption('outputNetcdfCaseName') , '.', num2str(ind),...
    '-', num2str( phys.Nmesh ), '.nc' ];

intervalType = phys.getOption('outputIntervalType');
switch intervalType
    case NdgIOIntervalType.DeltaTime
        interval = phys.getOption('outputTimeInterval');
    case NdgIOIntervalType.DeltaStep
        interval = phys.getOption('outputStepInterval');
    otherwise
        msgID = [mfilename, ':outputIntervalTypeInvalid'];
        msgtext = 'The output type is unknow.';
        throw( MException(msgID, msgtext) );
end
ncfile = getOutputFile( NdgIOFileType.NetCDF, intervalType, filename, ...
    [dimTime, dimK, dimNp, dimNfield], [varTime, varField], interval);
end% func