m=0.8;
// set points
Point(1) = {0,0,0,m};
Point(2) = {0,30,0,m};
Point(3) = {10,30,0,m};
Point(4) = {40,25.2753,0,m};
Point(5) = {40,0,0,m};

// connect boundary
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 1};
Line Loop(1) = {-1, -2, -3, -4, -5};
Plane Surface(1) = {1};

// set boundary condition
Physical Line(5) = {1};
Physical Line(4) = {4};
Physical Line(2) = {2,3,5};
Physical Surface(1) = {1};

// set No. of elements on each lines
npx = 50;
npy = Ceil(npx*0.75);
npx2 = Ceil(npx/4);
npx3 = npx - npx2;
// devided the line with uniform elements
Transfinite Line{1} = npy+1;
Transfinite Line{4} = npy+1;
Transfinite Line{5} = npx+1;
Transfinite Line{2} = npx2+1;
Transfinite Line{3} = npx3+1;
// Define the Surface as transfinite, by specifying the 
// four corners of the transfinite interpolation
Transfinite Surface{1} = {1, 5, 4, 2};
// generate quadrilateral elements
// Recombine Surface{1};
