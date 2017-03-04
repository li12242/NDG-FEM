tol = 1e-10;

%% Quad: N = 3
N = 3;
quad = StdRegions.Quad(N);
Dr_ext = load('DrDs_Test/Dr_3.cc');
Ds_ext = load('DrDs_Test/Ds_3.cc');

for i = 1:quad.nNode^2
    assert( abs(quad.Dr(i) - Dr_ext(i)) <= tol);
    assert( abs(quad.Ds(i) - Ds_ext(i)) <= tol);
end


%% Quad: N = 4
N = 4;
quad = StdRegions.Quad(N);
Dr_ext = load('DrDs_Test/Dr_4.cc');
Ds_ext = load('DrDs_Test/Ds_4.cc');

for i = 1:quad.nNode^2
    assert( abs(quad.Dr(i) - Dr_ext(i)) <= tol);
    assert( abs(quad.Ds(i) - Ds_ext(i)) <= tol);
end

%% Quad: N = 5
N = 5;
quad = StdRegions.Quad(N);
Dr_ext = load('DrDs_Test/Dr_5.cc');
Ds_ext = load('DrDs_Test/Ds_5.cc');

for i = 1:quad.nNode^2
    assert( abs(quad.Dr(i) - Dr_ext(i)) <= tol);
    assert( abs(quad.Ds(i) - Ds_ext(i)) <= tol);
end