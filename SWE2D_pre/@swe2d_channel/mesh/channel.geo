lland = 500;
lopen = 500;
Point(1) = {0, -250, 0, lopen};
Point(2) = {20000, -250, 0, lopen};
Point(3) = {20000, 250, 0, lopen};
Point(4) = {0, 250, 0, lopen};
Line(1) = {1, 2};
Line(2) = {2, 3}; // east boundary
Line(3) = {3, 4}; 
Line(4) = {4, 1}; // west boundary
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

nx = 200;
ny = 3;
Transfinite Line{1} = nx+1;
Transfinite Line{2} = ny+1;
Transfinite Line{3} = nx+1;
Transfinite Line{4} = ny+1;

Transfinite Surface{1} = {1,2,4,3};
Recombine Surface{1};
Mesh.Smoothing = 10;

Physical Line(2) = {1, 3};
Physical Line(5) = {2}; // east bc - clamped
Physical Line(4) = {4}; // west bc - zero gradient
Physical Surface(1) = {1};
