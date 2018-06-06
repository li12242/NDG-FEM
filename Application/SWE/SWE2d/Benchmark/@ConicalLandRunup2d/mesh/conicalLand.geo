x1 = 0; x2 = 25;
y1 = 0; y2 = 30;
lr = 1;

Point(1) = {x1, y1, 0, lr}; 
Point(2) = {x2, y1, 0, lr};
Point(3) = {x2, y2, 0, lr};
Point(4) = {x1, y2, 0, lr};

Line(1) = {1,2}; Line(2) = {2,3};
Line(3) = {3,4}; Line(4) = {4,1};

xc = 12.5;
yc = 15.0;
r1 = 1.1;
r2 = 3.6;
le = 0.2;

// inner circle
xc11 = xc - r1; yc11 = yc - r1;
xc12 = xc + r1; yc12 = yc + r1;

Point(5) = {xc, yc, 0, le};
Point(6) = {xc, yc11, 0, le};
Point(7) = {xc12, yc, 0, le};
Point(8) = {xc, yc12, 0, le};
Point(9) = {xc11, yc, 0, le};

Circle(5) = {6, 5, 7}; Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9}; Circle(8) = {9, 5, 6};

// outer circle
xc21 = xc - r2; yc21 = yc - r2;
xc22 = xc + r2; yc22 = yc + r2;

Point(10) = {xc, yc21, 0, le};
Point(11) = {xc22, yc, 0, le};
Point(12) = {xc, yc22, 0, le};
Point(13) = {xc21, yc, 0, le};

Circle(9) = {10, 5, 11}; Circle(10) = {11, 5, 12};
Circle(11) = {12, 5, 13}; Circle(12) = {13, 5, 10};


Line Loop(1) = {1,2,3,4}; 
Line Loop(2) = {5,6,7,8};
Line Loop(3) = {9,10,11,12};

Plane Surface(1) = {1, 3};
Plane Surface(2) = {2, 3};
Plane Surface(3) = {2};

Physical Line(6) = {4}; // 指定水位边界
Physical Line(4) = {1,2,3}; // 可滑移固边界
Physical Surface(1) = {1,2,3};

//Mesh.Smoothing = 20;
//Recombine Surface{1,2,3,4};
