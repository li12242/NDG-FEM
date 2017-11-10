function [ E, G ] = matEvaluateFlux( obj, mesh, fphys )

[ E, G ] = mxEvaluateFlux2d( obj.hmin, obj.gra, mesh.EToR, fphys );

end