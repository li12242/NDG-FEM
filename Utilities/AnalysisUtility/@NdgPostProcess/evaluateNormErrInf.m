function [ err ] = evaluateNormErrInf( obj, fphys, fext )

err = zeros(obj.Nvar, 1);
for m = 1:obj.Nmesh
    for fld = 1:obj.Nvar
        temp = fphys{m}(:,:,fld) - fext{m}(:,:,fld);
        maxErr = max( abs( temp(:) ) );
        err(fld) = max( err(fld), maxErr );
    end
end
end% func