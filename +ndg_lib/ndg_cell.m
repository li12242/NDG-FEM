function [ stdcell ] = ndg_cell( N, type )
%NDG_CELL 根据输出单元类型与阶数返回对应单元对象。
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
end

