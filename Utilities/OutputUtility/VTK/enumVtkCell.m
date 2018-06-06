classdef enumVtkCell < int8
    
    enumeration
        VTK_EMPTY_CELL          (0)
        VTK_VERTEX              (1)
        VTK_POLY_VERTEX         (2)
        VTK_LINE                (3)
        VTK_POLY_LINE           (4)
        VTK_TRIANGLE            (5)
        VTK_TRIANGLE_STRIP      (6)
        VTK_POLYGON             (7)
        VTK_PIXEL               (8)
        VTK_QUAD                (9)
        VTK_TETRA               (10)
        VTK_VOXEL               (11)
        VTK_HEXAHEDRON          (12)
        VTK_WEDGE               (13)
        VTK_PYRAMID             (14)
        VTK_PENTAGONAL_PRISM    (15)
        VTK_HEXAGONAL_PRISM     (16)
    end
    
end