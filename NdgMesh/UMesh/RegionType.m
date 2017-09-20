%> @brief Define the elemental region types.
%
%> The RegionType defines the properties of each element, 
% ======================================================================
%> This class is part of the NDG-FEM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef RegionType < int8
   
    enumeration
        %> normal element
        Normal (0)
        %> sponge elements
        Sponge (1)
        %> refined elements
        Refine (2)
        %> wet elements (SWE solver)
        Wet (3)
        %> dry elements (SWE solver)
        Dry (4)
    end
    
end

