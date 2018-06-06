function outputResult( obj, time, field )
    % output time
    startInd = obj.outputStep;
    countInd = 1;
    netcdf.putVar(obj.ncfile.ncid, obj.timeVarableId, startInd, countInd, time);
    
    % output physical field
    startInd = [ 0, 0, 0, obj.outputStep ];
    countInd = [ obj.mesh.cell.Np, obj.mesh.K, obj.Nfield, 1 ];
    netcdf.putVar(obj.ncfile.ncid, obj.fieldVarableId, startInd, countInd, field);

    % increase output step num
    obj.outputStep = obj.outputStep + 1;
end
