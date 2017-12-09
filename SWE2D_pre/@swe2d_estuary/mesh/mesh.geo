dx = 400e3;
dx1 = Ceil(dx/2);
dy = 200e3;
dy1 = Ceil(dy/5);

lc = 10e3;
lr = lc*2.5;
Point(1) = {0.0, -dy1,  0, lc}; 
Point(2) = {0.0,  dy1,  0, lc};
Point(3) = {dx1, -dy1, 0, lc};
Point(4) = {dx1,  dy1, 0, lc};
//Point(5) = {3.4,   1.701, 0, lr}; // Okushiri island
x0 = dx1; y1 = dy1*2; r = dy1;
Point(11) = {x0, -y1, 0, lc}; // ch5
Point(12) = {x0,  y1, 0, lc}; // ch7

Point(5) = {dx1+dy1, -2*dy1,  0, lc}; 
Point(6) = {dx1+dy1,  2*dy1,  0, lc};
Point(7) = {dx1+dy1, -dy,  0, lc};
Point(8) = {dx1+dy1,  dy,  0, lc};
Point(9) = {dx, -dy, 0, lr};
Point(10) = {dx, dy, 0, lr};

Line(1) = {1,2}; 
Line(2) = {2,4}; Line(3) = {1,3}; 
Line(4) = {2,4}; Line(5) = {1,3}; 

Line(6) = {6,8}; Line(7) = {5,7}; 
Line(8) = {8,10}; Line(9) = {7,9}; 
Line(10) = {9,10}; 

// // circle around ch5
Circle(12) = {4, 12, 6}; 
Circle(11) = {3, 11, 5};
Line(13) = {4,3};

Line Loop(1) = {3, -13, -2, -1};
Line Loop(2) = {13, 11, 7, 9, 10, -8, -6, -12};
Plane Surface(1) = {1};
Plane Surface(2) = {2};

Physical Line(2) = {1,2,3,4,5,6,7,8,9,10,11,12}; // 可滑移固边界
Physical Surface(0) = {1,2};

npx = Ceil(dx1/lc/1.2);
npy = Ceil(dy1/lc*3);
Transfinite Line{2} = npx+1;
Transfinite Line{3} = npx+1;
Transfinite Line{1} = npy+1;
Transfinite Line{13} = npy+1;
Transfinite Surface{1} = {1, 3, 4, 2};

Mesh.Smoothing = 20;
Recombine Surface{1,2};