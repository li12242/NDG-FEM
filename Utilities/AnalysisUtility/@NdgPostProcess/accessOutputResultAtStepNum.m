function [ fphys ] = accessOutputResultAtStepNum(obj, stepId)

fphys = cell( obj.Nmesh, 1 );
for m = 1:obj.Nmesh
    Np = obj.meshUnion(m).cell.Np;
    K = obj.meshUnion(m).K;
    fphys{m} = ncread( obj.outputFile{m}, 'fphys', [1, 1, 1, stepId], [Np, K, obj.Nvar, 1]);
end

end% func