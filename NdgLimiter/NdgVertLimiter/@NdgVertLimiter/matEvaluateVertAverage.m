function [ fvert, fvmin, fvmax, cvar ] = matEvaluateVertAverage( obj, fphys, fieldId )

Nv = obj.Nv;
Nmesh = obj.Nmesh;

% calculate the cell averages
cvar = cell( Nmesh, 1 );
for m = 1:Nmesh
	cvar{m} = obj.meshUnion(m).GetMeshAverageValue( fphys{m}(:,:,fieldId) );
    %cvar{m} = obj.meshUnion(m).cell_mean( fphys{m}(:,:,fieldId) );
end

[ fvert, fvmin, fvmax ] = mxEvaluateVertAverage( cvar, Nv, obj.Nvc, obj.VToM, obj.VToK, obj.VToW );

end% func