function initPhysFromOptions( obj, mesh )

% call the superclass methods
initPhysFromOptions@NdgPhysMat( obj, mesh );
obj.matUpdateWetDryState( obj.fphys );

obj.zGrad = cell( obj.Nmesh, 1 );
for m = 1:obj.Nmesh
    mesh = obj.meshUnion(m);
    obj.zGrad{m} = zeros( mesh.cell.Np, mesh.K, 2 );
    zr = mesh.cell.Dr * obj.fphys{m}(:,:,3);
    obj.zGrad{m}(:,:,1) = mesh.rx .* zr;
end

%Friction Term
if obj.option.isKey('FrictionType') % the option exist
    switch obj.getOption('FrictionType')
        case SWEFrictionType.None
            obj.frictionSolver = NonFrictionTermSolver();
        case SWEFrictionType.Linear
            t = obj.getOption('FrictionCoefficient_r');
            obj.frictionSolver = LinearFrictionTermSolver1d(t);
        case SWEFrictionType.Manning
            n = obj.getOption('FrictionCoefficient_n');
            obj.frictionSolver = ManningFrictionSolver1d(n);
    end
else % the option does not exist
    obj.frictionSolver = NonFrictionTermSolver();
end

% Numerical flux
if obj.option.isKey('NumFluxType')
    if obj.getOption('NumFluxType') == SWENumFluxType.HLL
        obj.numfluxSolver = SWEHLLNumFluxSolver1d( );
    elseif( obj.getOption('NumFluxType') == SWENumFluxType.LF )
        obj.numfluxSolver = SWELFNumFluxSolver1d( );
    elseif( obj.getOption('NumFluxType') == SWENumFluxType.ROE )
        obj.numfluxSolver = SWERoeNumFluxSolver1d( );
    end
else
    obj.numfluxSolver = SWEHLLNumFluxSolver1d( );
end

if obj.option.isKey('SWELimiterType')
    switch obj.getOption('SWELimiterType')
        case SWELimiterType.OnDepth
            obj.limiterSolver = SWEDepthLimiter1d();
        case SWELimiterType.OnElevation
            obj.limiterSolver = SWEElevationLimiter1d();
    end
else
    obj.limiterSolver = SWEDepthLimiter1d();
end

end% func



