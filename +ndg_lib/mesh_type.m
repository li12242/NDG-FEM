classdef mesh_type < int8
    %MESH_TYPE 计算单元枚举类型
    %   Detailed explanation goes here
    
    enumeration
        Normal      (0) % 普通单元
        Sponge      (1) % sponge cell
        Refine      (2) % 细分单元
        Coarse      (3) % 粗单元已细分
        Wet         (4) % 湿单元
        Dry         (5) % 干单元
    end
    
    methods
    end
    
end

