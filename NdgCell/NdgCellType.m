%> @brief enumeration for standard cell type.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDGOM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef NdgCellType < int8
    enumeration
        Point   (0)
        Line    (1)
        Tri     (2)
        Quad    (3)
        PrismTri (4)
        PrismQuad (5)
    end
end% classdef