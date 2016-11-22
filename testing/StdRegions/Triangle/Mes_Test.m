% test of Vandermonde matrix
tol = 1e-12;

%% Tri: N = 3
Mes_ext = load('Mes_Test/Mes_3.cc');

% triangle shape
N = 3;
tri = StdRegions.Triangle(N);

for i = 1:tri.nNode*tri.nFaceNode
    assert( abs(tri.Mes(i) - Mes_ext(i)) <= tol);
end

%% Tri: N = 4
Mes_ext = load('Mes_Test/Mes_4.cc');

% triangle shape
N = 4;
tri = StdRegions.Triangle(N);

for i = 1:tri.nNode*tri.nFaceNode
    assert( abs(tri.Mes(i) - Mes_ext(i)) <= tol);
end

%% Tri: N = 5
Mes_ext = load('Mes_Test/Mes_5.cc');

% triangle shape
N = 5;
tri = StdRegions.Triangle(N);

for i = 1:tri.nNode*tri.nFaceNode
    assert( abs(tri.Mes(i) - Mes_ext(i)) <= tol);
end