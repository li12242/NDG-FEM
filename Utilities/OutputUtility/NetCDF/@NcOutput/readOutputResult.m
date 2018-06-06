function [ field ] = readOutputResult( obj, step )

if ( step > obj.outputStep )
    error( ['The output step number is ',  num2str(obj.outputStep), ...
        ' less than input step number ', num2str(step) , '.\n'] )
end

Np = obj.mesh.cell.Np;
K = obj.mesh.K;
field = ncread( obj.filename, 'fphys', [1, 1, 1, step], [Np, K, obj.Nfield, 1]);

end% func