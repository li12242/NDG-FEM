function matEvaluateCellSurfFlux( obj, fphys, fext )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for m = 1:obj.Nmesh
    mesh = obj.meshUnion(m);
    
    f1 = fphys{m}( mesh.eidM );
    f2 = fphys{m}( mesh.eidP );
    flux = evaluateNumericalFlux( f1, f2, mesh.nx, mesh.ny, obj.u0, obj.v0 );
    obj.frhs{m} = obj.frhs{m} - mesh.cell.SLIFT*(mesh.Js .* flux);
    
end

end
