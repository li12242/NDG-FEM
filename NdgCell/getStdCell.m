%======================================================================
%> @brief warp function to get different standard cell objects
%>
%> @param N order of basis functions 
%> @param type standard cell type
%>
%> @retval stdCell standard cell class.
%======================================================================
%> This function is part of the NDGOM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
%======================================================================
function [ stdCell ] = getStdCell( N, type )
%get_std_cell std cell warp function to get the specific cell object
%   Detailed explanation goes here

switch type
    case NdgCellType.Point
        stdCell = StdPoint(N);
    case NdgCellType.Line
        stdCell = StdLine(N);
    case NdgCellType.Tri
        stdCell = StdTri(N);
    case NdgCellType.Quad
        stdCell = StdQuad(N);
end
end% func

