function [ err ] = evaluateNromErrInf( obj )
err = zeros(obj.Nfield, 1);
for m = 1:obj.Nmesh
    for fld = 1:obj.Nfield
        temp = obj.fphys{m}(:,:,fld) - obj.fext{m}(:,:,fld);
        maxErr = max( abs( temp(:) ) );
        err(fld) = max( err(fld), maxErr );
    end
end
end% func