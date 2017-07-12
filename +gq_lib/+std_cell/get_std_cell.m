function [ stdcell ] = get_std_cell( N, type )
%GET_STD_CELL Summary of this function goes here
%   Detailed explanation goes here

switch type
    case ndg_lib.std_cell_type.Point
        stdcell = gq_lib.std_cell.point(N);
    case ndg_lib.std_cell_type.Line
        stdcell = gq_lib.std_cell.line(N);
    case ndg_lib.std_cell_type.Tri
        stdcell = gq_lib.std_cell.tri(N);
    case ndg_lib.std_cell_type.Quad
        stdcell = gq_lib.std_cell.quad(N);
end
end

