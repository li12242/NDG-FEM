function [ fphys ] = matEvaluatePostFunc( obj, fphys )

for m = 1:obj.Nmesh
    hc = obj.meshUnion(m).GetMeshAverageValue( fphys{m}(:,:,1) );
    qxc = obj.meshUnion(m).GetMeshAverageValue( fphys{m}(:,:,2) );
    qyc = obj.meshUnion(m).GetMeshAverageValue( fphys{m}(:,:,3) );
    fphys{m}(:,:,1:3) = mxEvaluatePostFunc2d( obj.hmin, fphys{m}, hc, qxc, qyc );
end

obj.matUpdateWetDryState( fphys );

end
