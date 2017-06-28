c = 0.1;
Point(1) = {-1, -1, 0, c};
Point(2) = { 1, -1, 0, c};
Point(3) = { 1,  1, 0, c};
Point(4) = {-1,  1, 0, c};
Line(1) = {1, 2};
Line(2) = {2, 3}; // east boundary
Line(3) = {3, 4}; 
Line(4) = {4, 1}; // west boundary
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Recombine Surface{1};
Mesh.Smoothing = 10;

Physical Line(5) = {1,2,3,4}; // east bc - clamped
Physical Surface(1) = {1};Coherence;
