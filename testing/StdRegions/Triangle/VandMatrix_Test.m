% test of Vandermonde matrix
tol = 1e-12;

%% Tri: N = 3
V_ext = load('VandMatrix_Test/VandMatrix_3.cc');

% triangle shape
N = 3;
tri = StdRegions.Triangle(N);

for i = 1:tri.nNode^2
    assert( abs(tri.V(i) - V_ext(i)) <= tol);
end

%% Tri: N = 4
V_ext = load('VandMatrix_Test/VandMatrix_4.cc');

% triangle shape
N = 4;
tri = StdRegions.Triangle(N);

for i = 1:tri.nNode^2
    assert( abs(tri.V(i) - V_ext(i)) <= tol);
end

%% Tri: N = 5
V_ext = load('VandMatrix_Test/VandMatrix_5.cc');

% triangle shape
N = 5;
tri = StdRegions.Triangle(N);

for i = 1:tri.nNode^2
    assert( abs(tri.V(i) - V_ext(i)) <= tol);
end