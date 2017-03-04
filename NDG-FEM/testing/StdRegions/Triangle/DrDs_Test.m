tol = 1e-10;

%% Tri: N = 3
N = 3;
tri = StdRegions.Triangle(N);
Dr_ext = load('DrDs_Test/Dr_3.cc');
Ds_ext = load('DrDs_Test/Ds_3.cc');

for i = 1:tri.nNode^2
    assert( abs(tri.Dr(i) - Dr_ext(i)) <= tol);
    assert( abs(tri.Ds(i) - Ds_ext(i)) <= tol);
end


%% Tri: N = 4
N = 4;
tri = StdRegions.Triangle(N);
Dr_ext = load('DrDs_Test/Dr_4.cc');
Ds_ext = load('DrDs_Test/Ds_4.cc');

for i = 1:tri.nNode^2
    assert( abs(tri.Dr(i) - Dr_ext(i)) <= tol);
    assert( abs(tri.Ds(i) - Ds_ext(i)) <= tol);
end

%% Tri: N = 5
N = 5;
tri = StdRegions.Triangle(N);
Dr_ext = load('DrDs_Test/Dr_5.cc');
Ds_ext = load('DrDs_Test/Ds_5.cc');

for i = 1:tri.nNode^2
    assert( abs(tri.Dr(i) - Dr_ext(i)) <= tol);
    assert( abs(tri.Ds(i) - Ds_ext(i)) <= tol);
end