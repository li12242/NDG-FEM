lopen = 500;
len = 20e3; // length
width = 500; // width
Point(1) = {0, -width/2, 0, lopen};
Point(2) = {len, -width/2, 0, lopen};
Point(3) = {len, width/2, 0, lopen};
Point(4) = {0, width/2, 0, lopen};
spglen = 500;
Point(5) = {len - spglen, -width/2, 0, lopen };
Point(6) = {len - spglen, width/2, 0, lopen };
Line(1) = {1, 5};
Line(2) = {5, 2};
Line(3) = {2, 3}; // east boundary
Line(4) = {3, 6};
Line(5) = {6, 4}; 
Line(6) = {4, 1}; // west boundary
Line(7) = {5, 6};

Line Loop(1) = {1, 7, 5, 6};
Line Loop(2) = {-7, 2, 3, 4};
Plane Surface(0) = {1}; // normal region
Plane Surface(1) = {2}; // sponge region

nx = 100;
ny = 3;
Transfinite Line{1} = Ceil(nx - nx*spglen/len)+1;
Transfinite Line{7} = ny+1;
Transfinite Line{5} = Ceil(nx - nx*spglen/len)+1;
Transfinite Line{6} = ny+1;
Transfinite Surface{0} = {1,5,6,4};

Transfinite Line{2} = Ceil(nx*spglen/len)+1;
Transfinite Line{3} = ny+1;
Transfinite Line{4} = Ceil(nx*spglen/len)+1;
Transfinite Surface{1} = {5,2,6,3};

Recombine Surface{0,1};
Mesh.Smoothing = 10;

Physical Line(2) = {1, 2, 5, 4};
Physical Line(5) = {3}; // east bc - clamped
Physical Line(4) = {6}; // west bc - zero gradient
Physical Surface(0) = {0};
Physical Surface(1) = {1};
