% test of the LGL node set with exact coordinate
tol = 10e-10;

%% Test 1: N=3
r_ext = load('Coor_Test/r_3.cc');
s_ext = load('Coor_Test/s_3.cc');

% triangle shape
N = 3;
tri = StdRegions.Triangle(N);

for i = 1:tri.nNode
    assert( abs(tri.r(i) - r_ext(i)) <= tol)
    assert( abs(tri.s(i) - s_ext(i)) <= tol)
end

%% Test 2: N=4
r_ext = load('Coor_Test/r_4.cc');
s_ext = load('Coor_Test/s_4.cc');

% triangle shape
N = 4;
tri = StdRegions.Triangle(N);

for i = 1:tri.nNode
    assert( abs(tri.r(i) - r_ext(i)) <= tol)
    assert( abs(tri.s(i) - s_ext(i)) <= tol)
end

%% Test 3: N=5
r_ext = load('Coor_Test/r_5.cc');
s_ext = load('Coor_Test/s_5.cc');

% triangle shape
N = 5;
tri = StdRegions.Triangle(N);

for i = 1:tri.nNode
    assert( abs(tri.r(i) - r_ext(i)) <= tol)
    assert( abs(tri.s(i) - s_ext(i)) <= tol)
end
