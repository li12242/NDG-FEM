
lc = .12;
lr = lc/3;
Point(1) = {0.0,   0.0,   0, lc}; 
Point(2) = {5.488, 0.0,   0, lc};
Point(3) = {0.0,   3.402, 0, lc};
Point(4) = {5.488, 3.402, 0, lc};
//Point(5) = {3.4,   1.701, 0, lr}; // Okushiri island
x0 = 4.521;
y1 = 1.196; y2 = 1.696; y3 = 2.196;
r = 0.1;
Point(6) = {x0, y1, 0, lr}; // ch5
Point(7) = {x0, y2, 0, lr}; // ch7
Point(8) = {x0, y3, 0, lr}; // ch9
// points around ch5
Point(9) = {x0 - r, y1, 0, lr}; Point(10) = {x0, y1 - r, 0, lr};
Point(11) = {x0 + r, y1, 0, lr}; Point(12) = {x0, y1 + r, 0, lr};

// points around ch7
Point(13) = {x0 - r, y2, 0, lr}; Point(14) = {x0, y2 - r, 0, lr};
Point(15) = {x0 + r, y2, 0, lr}; Point(16) = {x0, y2 + r, 0, lr};

// points around ch9
Point(17) = {x0 - r, y3, 0, lr}; Point(18) = {x0, y3 - r, 0, lr};
Point(19) = {x0 + r, y3, 0, lr}; Point(20) = {x0, y3 + r, 0, lr};

Line(1) = {1,2}; Line(2) = {2,4};
Line(3) = {4,3}; Line(4) = {3,1};

// circle around ch5
Circle(5) = {9, 6, 10}; Circle(6) = {10, 6, 11};
Circle(7) = {11, 6, 12}; Circle(8) = {12, 6, 9};
// circle around ch7
Circle(9) = {13, 7, 14}; Circle(10) = {14, 7, 15};
Circle(11) = {15, 7, 16}; Circle(12) = {16, 7, 13};
// circle around ch9
Circle(13) = {17, 8, 18}; Circle(14) = {18, 8, 19};
Circle(15) = {19, 8, 20}; Circle(16) = {20, 8, 17};

Line Loop(1) = {1,2,3,4}; 
Line Loop(5) = {5,6,7,8};
Line Loop(7) = {9,10,11,12};
Line Loop(9) = {13,14,15,16};
Plane Surface(1) = {1,5,7,9};
Plane Surface(2) = {5};
Plane Surface(3) = {7};
Plane Surface(4) = {9};

// 设定边界物理属性
Physical Line(6) = {4}; // 指定水位边界
Physical Line(2) = {1,2,3}; // 可滑移固边界
Physical Surface(1) = {1,2,3,4};

// 设定在岛屿及实测点附近加密

Mesh.Smoothing = 20;
//Recombine Surface{1,2,3,4};

