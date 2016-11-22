tol = 1e-10;

%% Quad: N = 3
V_ext = load('VandMatrix_Test/VandMatrix_3.cc');

N = 3;
quad = StdRegions.Quad(N);

for i = 1:quad.nNode^2
    assert( abs(quad.V(i) - V_ext(i)) <= tol);
end


%% Quad: N = 4
V_ext = load('VandMatrix_Test/VandMatrix_4.cc');

N = 4;
quad = StdRegions.Quad(N);

for i = 1:quad.nNode^2
    assert( abs(quad.V(i) - V_ext(i)) <= tol);
end

%% Quad: N = 5
V_ext = load('VandMatrix_Test/VandMatrix_5.cc');

N = 5;
quad = StdRegions.Quad(N);

for i = 1:quad.nNode^2
    assert( abs(quad.V(i) - V_ext(i)) <= tol);
end