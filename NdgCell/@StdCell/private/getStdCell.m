function [ stdCell ] = getStdCell( N, type )
%getStdCell std cell warp function to get the specific cell object
%   Detailed explanation goes here

switch type
    case enumStdCell.Point
        stdCell = StdPoint(N);
    case enumStdCell.Line
        stdCell = StdLine(N);
    case enumStdCell.Tri
        stdCell = StdTri(N);
    case enumStdCell.Quad
        stdCell = StdQuad(N);
end
end% func

