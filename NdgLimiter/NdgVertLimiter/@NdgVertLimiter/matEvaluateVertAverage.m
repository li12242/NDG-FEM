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

% for n = 1:Nv
%     Nk = obj.Nvc(n);
% 	for m = 1:Nk
% 		meshId = obj.VToM(m, n);
% 		cellId = obj.VToK(m, n);
% 		w = obj.VToW(m, n);
% 		fvert(n) = fvert(n) + w*cvar{meshId}(cellId);
%         fvmin(n) = min( fvmin(n), cvar{meshId}(cellId) );
%         fvmax(n) = max( fvmax(n), cvar{meshId}(cellId) );
% 	end
% end

end% func