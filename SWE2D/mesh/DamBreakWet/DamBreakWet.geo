Point(1) = {0, -100, 0, 4.8};
Point(2) = {0,  100, 0, 4.8};
Point(3) = {1000,  100, 0, 4.8};
Point(4) = {1000, -100, 0, 4.8};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

nx = 100;
ny = 40;
Transfinite Line{1} = ny+1;
Transfinite Line{2} = nx+1;
Transfinite Line{3} = ny+1;
Transfinite Line{4} = nx+1;

Transfinite Surface{1} = {1,2,4,3};

Physical Line(3) = {2,4};
Physical Line(4) = {1,3};
Physical Surface(1) = {1};
