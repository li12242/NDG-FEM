lopen = 500;
len = 20e3; // length
width = 500; // width
Point(1) = {0, -width/2, 0, lopen};
Point(2) = {len, -width/2, 0, lopen};
Point(3) = {len, width/2, 0, lopen};
Point(4) = {0, width/2, 0, lopen};
Line(1) = {1, 2};
Line(2) = {2, 3}; // east boundary
Line(3) = {3, 4}; 
Line(4) = {4, 1}; // west boundary

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(0) = {1}; // normal region

nx = 100;
ny = 3;
Transfinite Line{1} = nx+1;
Transfinite Line{2} = ny+1;
Transfinite Line{3} = nx+1;
Transfinite Line{4} = ny+1;
Transfinite Surface{0} = {1,2,4,3};

Recombine Surface{0};
Mesh.Smoothing = 10;

Physical Line(2) = {1, 3};
Physical Line(5) = {2}; // east bc - clamped
Physical Line(4) = {4}; // west bc - zero gradient
Physical Surface(0) = {0};
