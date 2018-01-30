function [ fvmin, fvmax, cvar ] = matEvaluateVertAverage( ...
    obj, ...
    fphys, ...
    fieldId ...
    )

Nv = obj.Nv;
Nmesh = obj.Nmesh;

% calculate the cell averages
cvar = cell( Nmesh, 1 );
for m = 1:Nmesh
	cvar{m} = obj.meshUnion(m).GetMeshAverageValue(...
        fphys{m}(:,:,fieldId) );
end

[ fvmin, fvmax ] = mxEvaluateVertAverage( ...
    cvar, Nv, obj.Ncv, obj.VToM, obj.VToK );

end

