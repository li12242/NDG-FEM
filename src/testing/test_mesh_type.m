function test_mesh_type

shape = femlib.element_type('dim',1, 'vertice', 2, 'degree', 6, 'name', 'line');
EToV = [1,2; 2,3]; nodes = 3; elements = 2; name = 'Coor';

% new mesh
mesh = femlib.mesh_type('dim',1, 'EToV', EToV, ...
    'vertice', nodes, 'element', elements, 'shape', shape, 'name', name);
coord = [1,2];

mesh.setVerticeCoordMatrix(coord)
mesh

end
