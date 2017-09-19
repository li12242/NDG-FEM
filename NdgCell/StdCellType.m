%> @brief Define the enumeriations of the standard cell types
% 
classdef StdCellType < int8
    
    enumeration
        %> StdPoint type
        Point   (0)
        %> StdLine type
        Line    (1)
        %> StdTri type
        Tri     (2)
        %> StdQuad type
        Quad    (3)
        %> StdPrismTri type
        PrismTri (4)
        %> StdPrismQuad type
        PrismQuad (5)
    end
end