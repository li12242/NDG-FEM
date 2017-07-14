
tol = 10e-10;

%% Line: N=1
Mes_ext = load('Mes-test/Mes_1.cc');

% triangle shape
N = 1;
Lin = StdRegions.Line(N);
for i = 1:numel(Mes_ext)
    assert( abs(Lin.Mes(i) - Mes_ext(i)) <= tol)
end

%% Line: N=2
Mes_ext = load('Mes-test/Mes_2.cc');

% triangle shape
N = 2;
Lin = StdRegions.Line(N);
for i = 1:numel(Mes_ext)
assert( abs(Lin.Mes(i) - Mes_ext(i)) <= tol)
end

%% Line: N=3
Mes_ext = load('Mes-test/Mes_3.cc');
% triangle shape
N = 3;
Lin = StdRegions.Line(N);
for i = 1:numel(Mes_ext)
assert( abs(Lin.Mes(i) - Mes_ext(i)) <= tol)
end
%% Line: N=4
Mes_ext = load('Mes-test/Mes_4.cc');
% triangle shape
N = 4;
Lin = StdRegions.Line(N);
for i = 1:numel(Mes_ext)
assert( abs(Lin.Mes(i) - Mes_ext(i)) <= tol)
end
%% Line: N=5
Mes_ext = load('Mes-test/Mes_5.cc');
% triangle shape
N = 5;
Lin = StdRegions.Line(N);
for i = 1:numel(Mes_ext)
assert( abs(Lin.Mes(i) -Mes_ext(i)) <= tol)
end