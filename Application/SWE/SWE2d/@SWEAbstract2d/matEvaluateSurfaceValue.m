function [ fM, fP ] = matEvaluateSurfaceValue( obj, mesh, fphys, fext )

[ fM, fP ] = mxEvaluateSurfaceValue(  ...
    obj.hmin, obj.gra, ...
    mesh.eidM, mesh.eidP, ...
    mesh.nx, mesh.ny, ...
    mesh.eidtype, ...
    fphys, fext, ...
    mesh.EToE, mesh.EToR );

end% function