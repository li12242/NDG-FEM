classdef bc_type < int8
    %BC_TYPE 开边界枚举类型
    %   Detailed explanation goes here
    
    enumeration
        Inner           (0)
        SlipWall        (2)
        NonSlipWall     (3)
        ZeroGrad        (4)
        Clamped         (5)
        ClampedDepth    (6)
        ClampedVel      (7)
        Flather         (8)
    end
end

