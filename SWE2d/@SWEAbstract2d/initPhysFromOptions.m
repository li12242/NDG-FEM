function initPhysFromOptions( obj, mesh )

% call the superclass methods
initPhysFromOptions@NdgPhysMat( obj, mesh );
obj.matUpdateWetDryState( obj.fphys );

for m = 1:obj.Nmesh
    mesh = obj.meshUnion(m);
    obj.fphys{m}(:,:,5) = ...
        mesh.rx .* (mesh.cell.Dr * obj.fphys{m}(:,:,4)) + ...
        mesh.sx .* (mesh.cell.Ds * obj.fphys{m}(:,:,4));
    obj.fphys{m}(:,:,6) = ...
        mesh.ry .* (mesh.cell.Dr * obj.fphys{m}(:,:,4)) + ...
        mesh.sy .* (mesh.cell.Ds * obj.fphys{m}(:,:,4));
end

%Coriolis Term
if obj.option.isKey('CoriolisType') % the option exist
    switch obj.getOption('CoriolisType')
        case CoriolisType.None
            obj.coriolisSolver = NonCoriolisTermSolver();
        case CoriolisType.Beta
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
        case FrictionType.None
            obj.frictionSolver = NonFrictionTermSolver();
        case FrictionType.Linear
            t = obj.getOption('FrictionCoefficient_r');
            obj.frictionSolver = LinearFrictionTermSolver(t);
    end
else % the option does not exist
    obj.frictionSolver = NonFrictionTermSolver();
end

%Wind Term
if obj.option.isKey('WindType') % the option exist
    switch obj.getOption('WindType')
        case WindType.None
            obj.windSolver = NonWindTermSolver();
        case WindType.Stress
            q = obj.getOption('DensityofWater');
            obj.windSolver = StressWindTermSolver(q);
            %             obj.windSolver.rou = q;
        case WindType.UV
            o = obj.getOption('WindSterssCoefficient_cd');
            p = obj.getOption('DensityofAir');
            q = obj.getOption('DensityofWater');
            obj.windSolver = UVWindTermSolver(o, p, q);
    end
else % the option does not exist
    obj.windSolver = NonWindTermSolver();
end


end% func



