function [ mass ] = checkMassVolume( obj, varId )

time = ncread( obj.outputFile{1}, 'time' );
mass = zeros( obj.Nt, 1 );
for t = 1:obj.Nt
    [ fphys ] = obj.accessOutputResultAtStepNum(t);
    for m = 1:obj.Nmesh
        %mesh = obj.meshUnion(m);
        %temp = mesh.cell.V \ fphys{m}(:, :, varId);
        %mass(t) = mass(t) + sum( temp(1, :) );
        mass( t ) = mass( t ) + sum( obj.meshUnion(m).GetMeshIntegralValue( fphys{m}(:,:,varId) ) );
    end
end
temp = (mass - mass(1))./mass(1);

plot( time, temp );

end

