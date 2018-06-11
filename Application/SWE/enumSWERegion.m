classdef enumSWERegion < int8
    %ENUMSWEREGION Summary of this class goes here
    %   Detailed explanation goes here
    
    enumeration
        Sponge              (3) % sponge cell
        Wet                 (4) % well cell (SWE)
        Dry                 (5) % dry cell (SWE)
        PartialWet          (6)
        PartialWetFlood     (7)
        PartialWetDamBreak  (8)
    end
    
end

