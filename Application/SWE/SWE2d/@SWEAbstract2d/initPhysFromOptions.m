function initPhysFromOptions( obj, mesh )

% call the superclass methods
initPhysFromOptions@NdgPhysMat( obj, mesh );
obj.matUpdateWetDryState( obj.fphys );

obj.zGrad = cell( obj.Nmesh, 1 );
for m = 1:obj.Nmesh
    mesh = obj.meshUnion(m);
    
    obj.zGrad{m} = zeros( mesh.cell.Np, mesh.K, 2 );
    zr = mesh.cell.Dr * obj.fphys{m}(:,:,4);
    zs = mesh.cell.Ds * obj.fphys{m}(:,:,4);
    obj.zGrad{m}(:,:,1) = mesh.rx .* zr + mesh.sx .* zs;
    obj.zGrad{m}(:,:,2) = mesh.ry .* zr + mesh.sy .* zs;
end

%Coriolis Term
[ obj.coriolisSolver ] = initCoriolisSolver( obj );
%Friction Term
[ obj.frictionSolver ] = initFrictionSolver( obj );
%Wind Term
[ obj.windSolver ] = initWindSolver( obj );
[ obj.numfluxSolver ] = initNumFluxSolver( obj );
[ obj.limiterSolver ] = initLimiter( obj );

end% func



