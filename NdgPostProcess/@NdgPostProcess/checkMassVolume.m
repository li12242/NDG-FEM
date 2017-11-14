function [ mass ] = checkMassVolume( obj, varId )

time = ncread( obj.outputFile{1}, 'time' );
mass = zeros( obj.Nt, 1 );
for t = 1:obj.Nt
    [ fphys ] = obj.accessOutputResultAtStepNum(t);
    for m = 1:obj.Nmesh
        mass( t ) = mass( t ) + sum( obj.meshUnion(m).GetMeshIntegralValue( fphys{m}(:,:,varId) ) );
    end
end
plot( time, mass );

end

