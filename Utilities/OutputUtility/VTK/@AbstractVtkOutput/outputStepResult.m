function outputStepResult( obj, step, field )
%myFun - Description
%
% Syntax: output = myFun(input)
%
% Long description

filename = [ obj.casename, '/', obj.casename, '.', ...
num2str(step, '%04d'), '.vtk' ];
fp = fopen(filename, 'w');

fprintf(fp, '# vtk DataFile Version 2\n');
fprintf(fp, ['NDG-FEM ', obj.casename, ' Simulation\n']);
fprintf(fp, 'ASCII\n');
fprintf(fp, 'DATASET UNSTRUCTURED_GRID\n');

fprintf(fp, '\nPOINTS %d double\n', obj.Npoint);
fprintf(fp, '%12.20f  %12.20f  %12.20f\n', obj.Points);

fprintf(fp, '\nCELLS %d %d\n', obj.Ncell, obj.Ncon);
dataFormat = [num2str(obj.Np, '%d'), ' ', repmat('%d ', 1, obj.Np), '\n'];
fprintf(fp, dataFormat, obj.CellVertList);

fprintf(fp, '\nCELL_TYPES %d\n', obj.Ncell);
dataFormat = [repmat('%d ', 1, 12), '\n'];
fprintf(fp, dataFormat, obj.ctype);

fprintf(fp, '\n\nPOINT_DATA %d\n', obj.Npoint);

fprintf(fp, 'SCALARS field%d double 1\n', 1);
fprintf(fp, 'LOOKUP_TABLE default\n');
fprintf(fp, '%20.12f \n', field);

fclose(fp);

end