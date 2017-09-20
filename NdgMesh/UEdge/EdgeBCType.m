%> @brief Here we have a brief description of the class.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDG-FEM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef EdgeBCType < int8

    enumeration
        %> The edge is refined and will not be involved in computation
        Refined         (0)
        %> The 
        Inner           (1)
        SlipWall        (2)
        NonSlipWall     (3)
        ZeroGrad        (4)
        Clamped         (5)
        ClampedDepth    (6)
        ClampedVel      (7)
        Flather         (8)
    end
    
end

