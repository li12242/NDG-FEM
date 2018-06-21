function ConnectMeshUnion( obj, meshId, meshUnion )

obj.ind = meshId;
obj.EToM = ones( obj.cell.Nface, obj.K ) * meshId;
for m = 1 : numel( meshUnion )
    if m == meshId
        continue; % jump this cycle
    end

    mesh1 = meshUnion(m);
    [ obj.EToM, obj.EToE, obj.EToF ] ...
        = mxAssembleMeshConnection( obj.ind, m, ...
        obj.K, mesh1.K, obj.cell.Nface, mesh1.cell.Nface, ...
        obj.cell.FToV, mesh1.cell.FToV, ...
        obj.EToV, mesh1.EToV, obj.EToM, obj.EToE, obj.EToF);
end

end% func