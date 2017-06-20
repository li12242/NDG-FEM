function [ limiter ] = BJ( mesh )
%BJ Construct the 
%   Detailed explanation goes here

switch mesh.cell.type
    case ndg_lib.std_cell_type.Line
        limiter = ndg_utility.limiter.BJ.BJ_line(mesh);
    case ndg_lib.std_cell_type.Tri
    case ndg_lib.std_cell_type.Quad
end

end

