%> @brief Get the specific StdCell object.
%>
%> The user should input the cell order and type.
%>
%> @param N order 
%> @param type standard cell type from StdCellType
%>
%> @retval stdcell standard cell object.
%======================================================================
%> This class is part of the NDG-FEM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
%======================================================================
function [ stdcell ] = getStdCell( N, types )

if (numel(N) ~= numel(types))
    msgID = 'getStdCell:';
    msgtext = 'The input numbers of the order and types should equal';
    ME = MException(msgID, msgtext);
    throw(ME);
end

for n = 1:numel(N)
    stdcell(n) = createStdCell( N(n), types(n) );
end

end% func

function cell = createStdCell( N, type )

switch type
    case StdCellType.Point
        cell = StdPoint(N);
    case StdCellType.Line
        cell = StdLine(N);
    case StdCellType.Tri
        cell = StdTri(N);
    case StdCellType.Quad
        cell = StdQuad(N);
    otherwise
        
end
end

