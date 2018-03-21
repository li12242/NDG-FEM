%> @brief emumeration for mesh region types.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDGOM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef NdgRegionType < int8
    
    enumeration
        Normal      (1) % 
        Refine      (2) % refined cell, not participate in calculation
        Sponge      (3) % sponge cell
        Wet         (4) % well cell (SWE)
        Dry         (5) % dry cell (SWE)
        PartialWet  (6)
        PartialWetFlood    (7)
        PartialWetDamBreak (8)
    end
end