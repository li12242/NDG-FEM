function initPhysFromOptions( obj, mesh )
    % call the superclass methods
    initPhysFromOptions@SWEAbstractCB2d( obj, mesh );
    % set the physical field for the NdgPhysMat solver

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
%             obj.windSolver.cd = o;
%             obj.windSolver.rouair = p;
%             obj.windSolver.rou = q;
    end
else % the option does not exist
    obj.windSolver = NonWindTermSolver();
end

end% func



