classdef DiffusionAbstract2d < NdgPhysMat
    %DIFFUSIONABSTRACT2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties ( Constant )
        %> number of physical field [ C, tau1, tau2 ]
        Nfield = 3
        %> number of variable field
        Nvar = 1;
        %> index of variable in physical field
        varFieldIndex = 1;
    end
    
    methods
    end
    
end

