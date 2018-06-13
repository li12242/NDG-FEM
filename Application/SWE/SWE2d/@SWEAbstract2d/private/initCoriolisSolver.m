function [ coriolisSolver ] = initCoriolisSolver( obj )
%INITCORIOLISSOLVER Summary of this function goes here
%   Detailed explanation goes here

if obj.option.isKey('CoriolisType') % the option exist
    switch obj.getOption('CoriolisType')
        case enumSWECoriolis.None
            coriolisSolver = NonCoriolisTermSolver();
        case enumSWECoriolis.Beta
            m = obj.getOption('CoriolisParameter_f0');
            n = obj.getOption('CoriolisParameter_beta');
            coriolisSolver = BetaApproCoriolisTermSolver(m, n);
        case enumSWECoriolis.Latitude
            file = obj.getOption('LatitudeFilePath');
            coriolisSolver = LatitudeCoriolisTermSolver(obj, file);
    end
else % the option does not exist
    coriolisSolver = NonCoriolisTermSolver();
end

end

