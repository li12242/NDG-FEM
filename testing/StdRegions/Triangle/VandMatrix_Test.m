% test of Vandermonde matrix
tol = 1e-12;

%% Test 1: N = 3
V_ext = load('VandMatrix_3.cc');

% triangle shape
N = 3;
tri = StdRegions.Triangle(N);

for i = 1:tri.nNode^2
    assert( abs(tri.V(i) - V_ext(i)) <= tol);
end

%% Test 2: N = 4

%% Test 3: N = 5