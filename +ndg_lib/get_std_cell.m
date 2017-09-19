function [ stdcell ] = get_std_cell( N, type )
%get_std_cell std cell warp function to get the specific cell object
%   Detailed explanation goes here

switch type
    case ndg_lib.std_cell_type.Point
        stdcell = ndg_lib.std_cell.point(N);
    case ndg_lib.std_cell_type.Line
        stdcell = ndg_lib.std_cell.line(N);
    case ndg_lib.std_cell_type.Tri
        stdcell = ndg_lib.std_cell.tri(N);
    case ndg_lib.std_cell_type.Quad
        stdcell = ndg_lib.std_cell.quad(N);
end
end% func

