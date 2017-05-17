classdef mesh_type < int8
    %MESH_TYPE Summary of this class goes here
    %   Detailed explanation goes here
    
    enumeration
        Normal      (0)
        Sponge      (1) % sponge cell
        Refine      (2) % refine cell
        Coarse      (3) % coares cell to be refined
        Wet         (4)
        Dry         (5)
    end
    
    methods
    end
    
end

