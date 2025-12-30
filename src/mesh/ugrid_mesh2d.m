meshfile = 'Application/SWE/SWE2d/Benchmark/@ConicalLandRunup2d/mesh/conicalLand.msh';
mesh = read_gmsh2d(meshfile);

vert = ugrid_set(mesh.Nv, "vert");
cell = ugrid_set(mesh.Ntri, "cell");

map_c2v = ugrid_map(cell, 1, vert, 1, )