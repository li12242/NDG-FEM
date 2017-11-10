function [ fphys ] = accessOutputResultAtStepNum(obj, stepId)

fphys = cell( obj.Nmesh, 1 );
for m = 1:obj.Nmesh
    Np = obj.meshUnion(m).cell.Np;
    K = obj.meshUnion(m).K;
    filename = [obj.getOption('outputNetcdfCaseName'), '.', ...
        num2str(m), '-', num2str(obj.Nmesh), '.nc'];
    fphys{m} = ncread(filename, 'fphys', ...
        [1, 1, 1, stepId], [Np, K, obj.Nfield, 1]);
end

end% func