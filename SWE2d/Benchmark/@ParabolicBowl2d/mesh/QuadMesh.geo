c = 4000/20;
xmin = -4000; xmax = 4000; xmean = 0;
ymin = -4000; ymax = 4000; ymean = 0;

Point(1) = {xmin, ymin, 0, c};
Point(2) = {xmax, ymin, 0, c};
Point(3) = {xmax, ymax, 0, c};
Point(4) = {xmin, ymax, 0, c};
Point(5) = {xmean, ymin, 0, c};
Point(6) = {xmean, ymax, 0, c};

Line(1) = {1, 5};
Line(2) = {5, 6};
Line(3) = {6, 4};
Line(4) = {4, 1};
Line(5) = {5, 2};
Line(6) = {2, 3};
Line(7) = {3, 6};

Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5, 6, 7, -2};
Plane Surface(1) = {1};
Plane Surface(2) = {2};

Recombine Surface{1};
Recombine Surface{2};
Mesh.Smoothing = 8;

Physical Line(5) = {1,5,6,3,4,7};

Physical Surface(1) = {1};
Physical Surface(2) = {2};

