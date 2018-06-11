function [ windSolver ] = initWindSolver( obj )
%INITWINDSOLVER Summary of this function goes here
%   Detailed explanation goes here

if obj.option.isKey('WindType') % the option exist
    switch obj.getOption('WindType')
        case enumSWEWind.None
            windSolver = NonWindTermSolver();
        case enumSWEWind.Stress
            q = obj.getOption('DensityofWater');
            windSolver = StressWindTermSolver(q);
        case enumSWEWind.UV
            o = obj.getOption('WindSterssCoefficient_cd');
            p = obj.getOption('DensityofAir');
            q = obj.getOption('DensityofWater');
            windSolver = UVWindTermSolver(o, p, q);
    end
else % the option does not exist
    windSolver = NonWindTermSolver();
end

end
