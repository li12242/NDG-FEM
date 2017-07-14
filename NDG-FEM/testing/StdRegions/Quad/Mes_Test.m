% test of Vandermonde matrix
tol = 1e-12;

%% Quad: N = 3
Mes_ext = load('Mes_Test/Mes_3.cc');

% triangle shape
N = 3;
quad = StdRegions.Quad(N);

for i = 1:quad.nNode*quad.nFaceNode
    assert( abs(quad.Mes(i) - Mes_ext(i)) <= tol);
end

%% Quad: N = 4
Mes_ext = load('Mes_Test/Mes_4.cc');

% triangle shape
N = 4;
quad = StdRegions.Quad(N);

for i = 1:quad.nNode*quad.nFaceNode
    assert( abs(quad.Mes(i) - Mes_ext(i)) <= tol);
end

%% Quad: N = 5
Mes_ext = load('Mes_Test/Mes_5.cc');

% triangle shape
N = 5;
quad = StdRegions.Quad(N);

for i = 1:quad.nNode*quad.nFaceNode
    assert( abs(quad.Mes(i) - Mes_ext(i)) <= tol);
end