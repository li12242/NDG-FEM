%> @brief Enumeration of the standard cell types.
%> @note PrismQuad is not implemented (StdPrismQuad never completed);
%> the enum member is kept for value stability.
% ======================================================================
%> This class is part of the NDG-FEM software.
% ======================================================================
classdef enumStdCell < int8
    enumeration
        Point       (0)
        Line        (1)
        Tri         (2)
        Quad        (3)
        PrismTri    (4)
        PrismQuad   (5)
    end
end% classdef