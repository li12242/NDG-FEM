
lc = .12;
Point(1) = {0.0,   0.0,   0, lc}; 
Point(2) = {5.488, 0.0,   0, lc};
Point(3) = {0.0,   3.402, 0, lc};
Point(4) = {5.488, 3.402, 0, lc};
Point(5) = {3.4,   1.701, 0, lc}; // Okushiri island

Line(1) = {1,2}; Line(2) = {2,4};
Line(3) = {4,3}; Line(4) = {3,1};
Line Loop(1) = {1,2,3,4}; Plane Surface(1) = {1};

// 设定边界物理属性
Physical Line(6) = {4}; // 指定水位边界
Physical Line(2) = {1,2,3}; // 可滑移固边界
Physical Surface(1) = {1};

// 设定在岛屿附近加密
Field[1] = Attractor;
Field[1].NodesList = {5};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc / 4;
Field[2].LcMax = lc;
Field[2].DistMin = 0.1;
Field[2].DistMax = 0.6;

Field[3] = Min;
Field[3].FieldsList = {1, 2};
Background Field = 3;
Mesh.Smoothing = 20;
Recombine Surface{1};
