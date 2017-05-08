classdef bc_type < double
    %BC_TYPE 开边界枚举类型
    %   Detailed explanation goes here
    
    enumeration
        Inner           (0)
        SlipWall        (2)
        NonSlipWall     (3)
        ZeroGrad        (4)
        Clamped         (5)
    end
end

