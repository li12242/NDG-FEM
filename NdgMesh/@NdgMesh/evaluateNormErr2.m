function [ err ] = evaluateNormErr2( obj, fphys, fext )

Nvar = size( fphys, 3 );
err = zeros( Nvar, 1 );

totalArea = sum( obj.LAV );
for fld = 1 : Nvar
    temp = fphys(:,:,fld) - fext(:,:,fld);
    integralErr = sqrt( sum( obj.GetMeshIntegralValue( temp.^2 ) ) );
    err(fld) = integralErr/totalArea;
end

end