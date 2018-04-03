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
if obj.option.isKey('CoriolisType') % the option exist
    switch obj.getOption('CoriolisType')
        case SWECoriolisType.None
            obj.coriolisSolver = NonCoriolisTermSolver();
        case SWECoriolisType.Beta
            m = obj.getOption('CoriolisParameter_f0');
            n = obj.getOption('CoriolisParameter_beta');
            obj.coriolisSolver = BetaApproCoriolisTermSolver(m, n);
    end
else % the option does not exist
    obj.coriolisSolver = NonCoriolisTermSolver();
end

%Friction Term
if obj.option.isKey('FrictionType') % the option exist
    switch obj.getOption('FrictionType')
        case SWEFrictionType.None
            obj.frictionSolver = NonFrictionTermSolver();
        case SWEFrictionType.Linear
            t = obj.getOption('FrictionCoefficient_r');
            obj.frictionSolver = LinearFrictionTermSolver2d(t);
        case SWEFrictionType.Manning
            n = obj.getOption('FrictionCoefficient_n');
            obj.frictionSolver = ManningFrictionSolver2d( n );
    end
else % the option does not exist
    obj.frictionSolver = NonFrictionTermSolver();
end

%Wind Term
if obj.option.isKey('WindType') % the option exist
    switch obj.getOption('WindType')
        case SWEWindType.None
            obj.windSolver = NonWindTermSolver();
        case SWEWindType.Stress
            q = obj.getOption('DensityofWater');
            obj.windSolver = StressWindTermSolver(q);
            %             obj.windSolver.rou = q;
        case SWEWindType.UV
            o = obj.getOption('WindSterssCoefficient_cd');
            p = obj.getOption('DensityofAir');
            q = obj.getOption('DensityofWater');
            obj.windSolver = UVWindTermSolver(o, p, q);
    end
else % the option does not exist
    obj.windSolver = NonWindTermSolver();
end

if obj.option.isKey('NumFluxType')
    if obj.getOption('NumFluxType') == SWENumFluxType.HLL
        obj.numfluxSolver = SWEHLLNumFluxSolver2d( );
    elseif( obj.getOption('NumFluxType') == SWENumFluxType.LF )
        obj.numfluxSolver = SWELFNumFluxSolver2d(  );
    end
else
    obj.numfluxSolver = SWEHLLNumFluxSolver2d( );
end

if obj.option.isKey('SWELimiterType')
    if obj.getOption('SWELimiterType') == SWELimiterType.OnDepth
        obj.limiterSolver = SWEDepthLimiter();
    elseif obj.getOption('SWELimiterType') == SWELimiterType.OnElevation
        obj.limiterSolver = SWEElevationLimiter();
    end
else
    obj.limiterSolver = SWEDepthLimiter();
end
end% func



