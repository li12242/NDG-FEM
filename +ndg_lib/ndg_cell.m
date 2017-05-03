function [ stdcell ] = ndg_cell( N, type )
%NDG_CELL Summary of this function goes here
%   Detailed explanation goes here

switch type
    case ndg_lib.std_cell_type.Point
        stdcell = ndg_lib.std_cell.std_point(N);
    case ndg_lib.std_cell_type.Line
        stdcell = ndg_lib.std_cell.std_line(N);
    case ndg_lib.std_cell_type.Tri
        stdcell = ndg_lib.std_cell.std_tri(N);
    case ndg_lib.std_cell_type.Quad
        stdcell = ndg_lib.std_cell.std_quad(N);
end
end

