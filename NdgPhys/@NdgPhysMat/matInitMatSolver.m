function matInitMatSolver( obj )


% setup the output file object
obj.outputNcFile = [];
for n = 1:obj.Nmesh
    [ ncfile ] = makeBasicOutputNetcdfFile( obj, n, obj.meshUnion(n) );
    switch obj.getOption('outputIntervalType')
        case NdgIntervalType.Constant
            obj.outputTimeInterval = obj.getOption('outputTimeInterval');
        case NdgIntervalType.DeltaTime
            obj.outputTimeInterval = obj.getOption('outputTimeInterval');
        case NdgIntervalType.DeltaStep
            obj.outputStepInterval = obj.getOption('outputStepInterval');
    end
    obj.outputNcFile = [obj.outputNcFile, ncfile];
end

%> set the final time
obj.ftime = obj.getOption('finalTime');
obj.outputStep = 0;

%> choose the limiter 
switch obj.getOption('LimiterType')
    case NdgLimiterType.None
    case NdgLimiterType.VertexLimiter
        obj.limiter = NdgVertLimiter2d(obj.meshUnion);
end
end% func

function [ ncfile ] = makeBasicOutputNetcdfFile( phys, ind, mesh )
dimTime = NdgNcDim('Nt', 0);
dimK = NdgNcDim('K', mesh.K);
dimNp = NdgNcDim('Np', mesh.cell.Np);
dimNfield = NdgNcDim('Nfield', phys.Nfield);

varTime = NdgNcVar('time', dimTime, NdgNcType.NC_DOUBLE );
varField = NdgNcVar('fphys', [dimNp, dimK, dimNfield, dimTime], NdgNcType.NC_DOUBLE);

filename = [phys.getOption('outputNetcdfCaseName') , '.', num2str(ind),...
    '-', num2str( phys.Nmesh ), '.nc'];
ncfile = NdgNcOutputFile(filename, [dimTime, dimK, dimNp, dimNfield], [varTime, varField]);
ncfile.defineIntoNetcdfFile(); % define the NetCDF file;
end% func