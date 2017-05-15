function [ limiter ] = BJ( mesh, cell )
%BJ Construct the 
%   Detailed explanation goes here

switch cell.type
    case ndg_lib.std_cell_type.Line
        limiter = ndg_utility.limiter.BJ.BJ_line(mesh, cell);
    case ndg_lib.std_cell_type.Tri
    case ndg_lib.std_cell_type.Quad
end

end

