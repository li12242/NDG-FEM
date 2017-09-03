m=0.8;
theta = Pi/6;
x1 = 0; x2 = 1e3;
y1 = 0; y2 = 20;
xc = x2/2; yc = y2/2;

cost = Cos(theta);
sint = Sin(theta);
// set points
Point(1) = {cost*x1 + sint*y1, -sint*x1 + cost*y1, 0,m};
Point(2) = {cost*x2 + sint*y1, -sint*x2 + cost*y1, 0,m};
Point(3) = {cost*x2 + sint*y2, -sint*x2 + cost*y2, 0,m};
Point(4) = {cost*x1 + sint*y2, -sint*x1 + cost*y2, 0,m};
Point(5) = {cost*x2/2 + sint*y1, -sint*x2/2 + cost*y1, 0,m};
Point(6) = {cost*x2/2 + sint*y2, -sint*x2/2 + cost*y2, 0,m};

// connect boundary
Line(1) = {1, 5};
Line(2) = {5, 2};
Line(3) = {2, 3};
Line(4) = {3, 6};
Line(5) = {6, 4};
Line(6) = {4, 1};
Line(7) = {5, 6};

Line Loop(1) = {7, 5, 6, 1};
Line Loop(2) = {7, -4, -3, -2};
Plane Surface(1) = {1};
Plane Surface(2) = {2};

// set boundary condition
Physical Line(4) = {3,6};
Physical Line(2) = {1,2,4,5};
Physical Surface(0) = {1};
Physical Surface(1) = {2};

// set No. of elements on each lines
npx = 300;
npy = 1;

// devided the line with uniform elements
Transfinite Line{1} = npx/2+1;
Transfinite Line{2} = npx/2+1;
Transfinite Line{3} = npy+1;
Transfinite Line{4} = npx/2+1;
Transfinite Line{5} = npx/2+1;
Transfinite Line{6} = npy+1;
Transfinite Line{7} = npy+1;

// Define the Surface as transfinite, by specifying the 
// four corners of the transfinite interpolation
Transfinite Surface{1} = {5, 6, 4, 1};
Transfinite Surface{2} = {5, 6, 2, 3};
// generate quadrilateral elements
Recombine Surface{1, 2};
