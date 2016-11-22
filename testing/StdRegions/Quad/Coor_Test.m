tol = 1e-10;

%% Quad: N = 3
r_ext = load('Coor_Test/r_3.cc');
s_ext = load('Coor_Test/s_3.cc');

N = 3;
quad = StdRegions.Quad(N);

for i = 1:quad.nNode
    assert( abs(quad.r(i) - r_ext(i)) <= tol)
    assert( abs(quad.s(i) - s_ext(i)) <= tol)
end

%% Quad: N = 4
r_ext = load('Coor_Test/r_4.cc');
s_ext = load('Coor_Test/s_4.cc');

N = 4;
quad = StdRegions.Quad(N);

for i = 1:quad.nNode
    assert( abs(quad.r(i) - r_ext(i)) <= tol)
    assert( abs(quad.s(i) - s_ext(i)) <= tol)
end

%% Quad: N = 5
r_ext = load('Coor_Test/r_5.cc');
s_ext = load('Coor_Test/s_5.cc');

N = 5;
quad = StdRegions.Quad(N);

for i = 1:quad.nNode
    assert( abs(quad.r(i) - r_ext(i)) <= tol)
    assert( abs(quad.s(i) - s_ext(i)) <= tol)
end