function initPhysFromOptions( obj, mesh )
    % call the superclass methods
    initPhysFromOptions@NdgPhysMat( obj, mesh );
    
    % set the bottom gradient
    for m = 1:obj.Nmesh
        mesh = obj.meshUnion(m);
        [ rx, ry, ~, sx, sy, ~, ~, ~, ~, ~ ] ...
                = mesh.cell.assembleJacobianMatrix( mesh.x, mesh.y, mesh.z );
        obj.fphys{m}(:,:,5) = rx .* (mesh.cell.Dr * obj.fphys{m}(:,:,4)) + ...
            sx .* (mesh.cell.Ds * obj.fphys{m}(:,:,4));
        obj.fphys{m}(:,:,6) = ry .* (mesh.cell.Dr * obj.fphys{m}(:,:,4)) + ...
            sy .* (mesh.cell.Ds * obj.fphys{m}(:,:,4));
        
        if ~obj.option.isKey('WellBlancedType')
            obj.fphys{m}(:,:,4) = 0; % set the bottom elevation z to 0
        end
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
    %             obj.coriolisSolver.f0 = m;
    %             obj.coriolisSolver.beta = n;
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
    %             obj.frictionSolver.r = t;
        end
    else % the option does not exist
        obj.frictionSolver = NonFrictionTermSolver();
    end

end% func



