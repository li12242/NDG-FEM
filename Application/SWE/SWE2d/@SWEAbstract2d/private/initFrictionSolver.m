function [ frictionSolver ] = initFrictionSolver( obj )
%INITFRICTIONSOLVER Summary of this function goes here
%   Detailed explanation goes here

if obj.option.isKey('FrictionType') % the option exist
    switch obj.getOption('FrictionType')
        case enumSWEFriction.None
            frictionSolver = NonFrictionTermSolver();
        case enumSWEFriction.Linear
            t = obj.getOption('FrictionCoefficient_r');
            frictionSolver = LinearFrictionTermSolver2d(t);
        case enumSWEFriction.Manning
            n = obj.getOption('FrictionCoefficient_n');
            frictionSolver = ManningFrictionSolver2d( n );
        case enumSWEFriction.Quadric
            t = obj.getOption('FrictionCoefficient_n');
            frictionSolver = QuadricFrictionTermSolver2d(t);
    end
else % the option does not exist
    frictionSolver = NonFrictionTermSolver();
end

end

